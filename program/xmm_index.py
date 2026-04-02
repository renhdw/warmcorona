#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import re
import subprocess
import tempfile
from datetime import datetime
import urllib.parse as up

import requests
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u


# =========================
# 0) 你的环境路径（按需改）
# =========================
XMM_BASE_DIR = Path("/home/hdw/data/XMM").resolve()
DEFAULT_OUT_ROOT = Path("~/data/monk/plot/warmcorona/test/2026.03.26/index").expanduser().resolve()

MODEL_SLAB = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_pure_Tom_slab_107.mod"
MODEL_SPHERE = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-1.0_5-25_smoothed_tv_Tom_new_107.mod"

WA_MRK509 = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/wa_Mrk509_xi_NH.mod"
WA_GENERIC = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/wa_cloudy_agn_v2wa_agn_xi_NH.mod"

XCM_SLAB_TEMPLATE = "fit_model_Tom_slab_{source}_{obsid}.xcm"
XCM_SPHERE_TEMPLATE = "fit_model_tv_Tom_sphere_{source}_{obsid}.xcm"

EXTRA_CMDS = [
    "plot ldata emo re",
    "freeze 3",
    "freeze 5",
    "new 8=2",
    "new 12=2",
    "new 14=4",
    "new 18=2",
    "freeze 15",
    "freeze 16",
    "freeze 17",
    "new 0",
]


# =========================
# 1) 正则
# =========================
_OBSID_RE = re.compile(r"^\d{10}$")
_CHI_XCM_RE = re.compile(r"^chi_\d+.*\.xcm$", re.IGNORECASE)

# 排除：chi 名字里含 gdiskwien 的都不要
_EXCLUDE_CHI_SUBSTRS = ["gdiskwien"]

# float
_FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
_SCI_RE = re.compile(r"([0-9.+-]+E[+-]?[0-9]+)", re.IGNORECASE)

# fit summary
_FITSTAT_RE = re.compile(r"Fit statistic\s*:\s*Chi-Squared\s+([0-9.]+)\s+using\s+(\d+)\s+bins", re.I)
_DOF_RE = re.compile(r"Null hypothesis probability of\s+([0-9.eE+-]+)\s+with\s+(\d+)\s+degrees of freedom", re.I)


# =========================
# 2) 数据结构
# =========================
@dataclass(frozen=True)
class ObsEntry:
    source: str
    obsid: str
    path: Path


@dataclass(frozen=True)
class XcmPick:
    path: Path
    mtime: float
    kind: str  # slab/sphere/unknown

    @property
    def mtime_str(self) -> str:
        return datetime.fromtimestamp(self.mtime).strftime("%Y-%m-%d %H:%M:%S")


# =========================
# 3) 小工具
# =========================
def _extract_sci(line: str) -> Optional[float]:
    m = _SCI_RE.search(line)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def _is_source_dir(p: Path) -> bool:
    return p.is_dir() and not p.name.startswith(".")


def _is_obsid_dir(p: Path) -> bool:
    return p.is_dir() and (_OBSID_RE.match(p.name) is not None)


def guess_chi_kind(name: str) -> str:
    n = name.lower()
    if "sphere" in n:
        return "sphere"
    if "slab" in n:
        return "slab"
    return "unknown"


def is_mrk509(source: str) -> bool:
    s = source.strip().lower().replace(" ", "")
    return s in {"mrk509", "markarian509"}


def patch_text_for_mrk509(source: str, text: str) -> str:
    if not is_mrk509(source):
        return text
    if WA_GENERIC in text and WA_MRK509 not in text:
        return text.replace(WA_GENERIC, WA_MRK509)
    return text


def should_exclude_chi(path: Path) -> bool:
    low = path.name.lower()
    return any(sub in low for sub in _EXCLUDE_CHI_SUBSTRS)


def compress_ranges(nums: List[int]) -> str:
    """把 [1,2,3,5,6,9] -> '1-3 5-6 9' (给 XSPEC error 用)"""
    if not nums:
        return ""
    nums = sorted(set(nums))
    parts = []
    start = prev = nums[0]
    for x in nums[1:]:
        if x == prev + 1:
            prev = x
            continue
        parts.append(f"{start}" if start == prev else f"{start}-{prev}")
        start = prev = x
    parts.append(f"{start}" if start == prev else f"{start}-{prev}")
    return " ".join(parts)


def parse_fit_summary(raw: str) -> Dict[str, object]:
    out = {"chi2": None, "bins": None, "nhp": None, "dof": None}
    m = _FITSTAT_RE.search(raw)
    if m:
        try:
            out["chi2"] = float(m.group(1))
            out["bins"] = int(m.group(2))
        except Exception:
            pass
    m2 = _DOF_RE.search(raw)
    if m2:
        try:
            out["nhp"] = float(m2.group(1))
            out["dof"] = int(m2.group(2))
        except Exception:
            pass
    return out


# =========================
# 4) 功能 1：扫描
# =========================
def scan_all_entries() -> List[ObsEntry]:
    out: List[ObsEntry] = []
    if not XMM_BASE_DIR.exists():
        return out
    for src_dir in sorted(XMM_BASE_DIR.iterdir(), key=lambda x: x.name.lower()):
        if not _is_source_dir(src_dir):
            continue
        src = src_dir.name
        for obs_dir in sorted(src_dir.iterdir(), key=lambda x: x.name):
            if not _is_obsid_dir(obs_dir):
                continue
            out.append(ObsEntry(source=src, obsid=obs_dir.name, path=obs_dir))
    return out


# =========================
# 5) 功能 2：挑 chi_xcm（优先 1 sphere + 1 slab；排除 gdiskwien）
# =========================
def list_chi_xcms(obs_dir: Path, recursive: bool) -> List[Path]:
    bases = [obs_dir, obs_dir / "xcm"]
    out: List[Path] = []
    for b in bases:
        if not b.exists() or not b.is_dir():
            continue
        if recursive:
            for p in b.rglob("*.xcm"):
                if p.is_file() and _CHI_XCM_RE.match(p.name) and (not should_exclude_chi(p)):
                    out.append(p)
        else:
            for p in b.iterdir():
                if p.is_file() and _CHI_XCM_RE.match(p.name) and (not should_exclude_chi(p)):
                    out.append(p)
    return out


def pick_chi_prefer_slab_sphere(obs_dir: Path, recursive: bool) -> Tuple[List[XcmPick], str]:
    """
    返回 picks(<=2) + note
    - 尽量返回：最新 sphere + 最新 slab
    - 如果只有一种：就返回该种最新两条（并在 note 里说明“不满足一球一板”）
    """
    allp: List[XcmPick] = []
    for f in list_chi_xcms(obs_dir, recursive):
        try:
            st = f.stat()
        except OSError:
            continue
        allp.append(XcmPick(path=f, mtime=st.st_mtime, kind=guess_chi_kind(f.name)))

    if not allp:
        return [], "NO chi_*.xcm found"

    slab = [p for p in allp if p.kind == "slab"]
    sph = [p for p in allp if p.kind == "sphere"]
    unk = [p for p in allp if p.kind == "unknown"]

    slab.sort(key=lambda x: x.mtime, reverse=True)
    sph.sort(key=lambda x: x.mtime, reverse=True)
    unk.sort(key=lambda x: x.mtime, reverse=True)

    if slab and sph:
        picks = [sph[0], slab[0]]
        picks.sort(key=lambda x: x.mtime, reverse=True)
        return picks, "OK (picked 1 sphere + 1 slab)"

    if slab:
        pool = slab + unk
        pool.sort(key=lambda x: x.mtime, reverse=True)
        return pool[:2], "WARN (only slab found; not 1 sphere + 1 slab)"
    if sph:
        pool = sph + unk
        pool.sort(key=lambda x: x.mtime, reverse=True)
        return pool[:2], "WARN (only sphere found; not 1 sphere + 1 slab)"

    unk.sort(key=lambda x: x.mtime, reverse=True)
    return unk[:2], "WARN (only unknown kind; not 1 sphere + 1 slab)"


def build_full_index(entries: List[ObsEntry], recursive: bool) -> Dict[Tuple[str, str], Tuple[List[XcmPick], str]]:
    idx: Dict[Tuple[str, str], Tuple[List[XcmPick], str]] = {}
    for e in entries:
        idx[(e.source, e.obsid)] = pick_chi_prefer_slab_sphere(e.path, recursive)
    return idx


def write_full_index_txt(
    out_root: Path,
    entries: List[ObsEntry],
    idx: Dict[Tuple[str, str], Tuple[List[XcmPick], str]],
) -> Path:
    out_dir = out_root / "xmm_index"
    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_txt = out_dir / f"xmm_full_index_{stamp}.txt"

    bysrc: Dict[str, List[ObsEntry]] = {}
    for e in entries:
        bysrc.setdefault(e.source, []).append(e)
    for s in bysrc:
        bysrc[s].sort(key=lambda x: x.obsid)
    bysrc = dict(sorted(bysrc.items(), key=lambda kv: kv[0].lower()))

    lines: List[str] = []
    lines.append("XMM FULL INDEX (always generated)")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"XMM_BASE_DIR: {XMM_BASE_DIR}")
    lines.append("Rule: prefer picking 1 sphere + 1 slab; exclude chi_* containing: " + ",".join(_EXCLUDE_CHI_SUBSTRS))
    lines.append("")

    for src, obslist in bysrc.items():
        lines.append(f"== {src} (obsids={len(obslist)}) ==")
        for e in obslist:
            picks, note = idx.get((e.source, e.obsid), ([], ""))
            lines.append(f"{src}/{e.obsid}  picked={len(picks)}  note={note}")
            if not picks:
                lines.append("  (NO chi_*.xcm found)")
            else:
                for k, p in enumerate(sorted(picks, key=lambda x: x.mtime, reverse=True), 1):
                    lines.append(f"  {k}) kind={p.kind:<7s}  {p.path}  mtime={p.mtime_str}")
        lines.append("")

    out_txt.write_text("\n".join(lines), encoding="utf-8")
    return out_txt


# =========================
# 6) 在线查 z/NH（Simbad + w3nh）
# =========================
def generate_name_candidates(name: str) -> List[str]:
    cand = [name.strip()]
    upper = name.strip().upper().replace(" ", "")
    if upper.startswith("IRASF"):
        rest = name.strip()[5:]
        cand += [f"IRAS F{rest}", f"IRAS {rest}"]
    elif upper.startswith("IRAS") and " " not in name.strip()[:5]:
        rest = name.strip()[4:]
        cand += [f"IRAS {rest}"]
    seen, out = set(), []
    for c in cand:
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out


def _simbad_row_to_radec_deg(row) -> Tuple[float, float]:
    ra_s = None
    dec_s = None
    for key in ("RA", "ra"):
        if key in row.colnames:
            ra_s = row[key]
            break
    for key in ("DEC", "dec"):
        if key in row.colnames:
            dec_s = row[key]
            break

    if hasattr(ra_s, "decode"):
        ra_s = ra_s.decode("utf-8")
    if hasattr(dec_s, "decode"):
        dec_s = dec_s.decode("utf-8")

    if isinstance(ra_s, str) and isinstance(dec_s, str):
        c = SkyCoord(ra_s, dec_s, unit=(u.hourangle, u.deg), frame="icrs")
        return float(c.ra.deg), float(c.dec.deg)

    try:
        return float(ra_s), float(dec_s)
    except Exception:
        raise RuntimeError(f"Cannot parse RA/DEC from SIMBAD row: RA={ra_s} DEC={dec_s}")


def query_simbad_basic(name: str) -> Dict[str, object]:
    custom = Simbad()
    custom.add_votable_fields("rvz_redshift")

    last_err = None
    result = None
    used_name = None
    for cand in generate_name_candidates(name):
        try:
            result = custom.query_object(cand)
        except Exception as e:
            last_err = e
            result = None
        if result is not None and len(result) > 0:
            used_name = cand
            break

    if result is None or len(result) == 0:
        raise RuntimeError(f"SIMBAD not found: {name}, last_err={last_err}")

    row = result[0]
    ra_deg, dec_deg = _simbad_row_to_radec_deg(row)

    z_val = None
    try:
        if "RVZ_REDSHIFT" in row.colnames:
            zv = row["RVZ_REDSHIFT"]
        else:
            zv = row["rvz_redshift"]
        if zv is not None:
            z_val = float(zv)
    except Exception:
        z_val = None

    main_id = name
    for key in ("MAIN_ID", "main_id"):
        if key in row.colnames:
            mid = row[key]
            if hasattr(mid, "decode"):
                mid = mid.decode("utf-8")
            main_id = str(mid)
            break

    return {"main_id": main_id, "ra_deg": ra_deg, "dec_deg": dec_deg, "z": z_val, "used_name": used_name or name}


def build_w3nh_params(ra_deg: float, dec_deg: float, radius_deg: float = 0.1) -> dict:
    entry = f"{ra_deg:.5f}, {dec_deg:.5f}"
    return {
        "Entry": entry,
        "NR": "GRB/SIMBAD+Sesame/NED",
        "CoordSys": "Equatorial",
        "equinox": "2000",
        "radius": f"{radius_deg}",
        "usemap": "0",
    }


def build_w3nh_url(params: dict) -> str:
    base = "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl"
    return f"{base}?{up.urlencode(params)}"


def fetch_w3nh_nh(ra_deg: float, dec_deg: float, radius_deg: float = 0.1) -> Dict[str, object]:
    base = "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl"
    params = build_w3nh_params(ra_deg, dec_deg, radius_deg)
    url = build_w3nh_url(params)

    headers = {"User-Agent": "Mozilla/5.0"}
    resp = requests.get(base, params=params, timeout=25, headers=headers)
    resp.raise_for_status()
    text = resp.text

    avg_nh, wavg_nh = None, None
    for line in text.splitlines():
        low = line.lower()
        if "average nh" in low and "weighted" not in low:
            v = _extract_sci(line)
            if v is not None:
                avg_nh = v
        if "weighted average nh" in low:
            v = _extract_sci(line)
            if v is not None:
                wavg_nh = v

    return {"avg_nh": avg_nh, "wavg_nh": wavg_nh, "url": url}


# =========================
# 7) 写 Tom xcm
# =========================
def build_params_for_source(z: Optional[float], nh_1e22: Optional[float]) -> List[float]:
    if z is None:
        z = 0.0
    if nh_1e22 is None:
        nh_1e22 = 0.01
    return [
        1.0, z, nh_1e22,      # zTBabs(nH,z), TBabs(nH)
        1.5, 100.0, 3e-3, 0,  # nthComp
        z, 1.0,               # nthComp z, norm
        0.5, 12.0, z, 1.0,    # warmcom te,tau,z,norm
        1.5, 1.0, 300.0, 0.0, z, 30.0, -1.0, 1.0  # xillver
    ]


def write_tom_xcm(obs_dir: Path, source: str, obsid: str, slab: bool, params: List[float], include_cpd: bool) -> Path:
    xcm_dir = obs_dir / "xcm"
    xcm_dir.mkdir(parents=True, exist_ok=True)

    model_path = MODEL_SLAB if slab else MODEL_SPHERE
    out_name = (XCM_SLAB_TEMPLATE if slab else XCM_SPHERE_TEMPLATE).format(source=source, obsid=obsid)
    out_path = xcm_dir / out_name

    param_lines = "\n".join(str(p) for p in params)
    extra_cmds = "\n".join(EXTRA_CMDS)
    cpd_line = "cpd /xs\n" if include_cpd else ""

    content = f"""
data PN_spectrum_grp.fits
{cpd_line}setpl energy
setpl add
ignore 1.8-2.2
ignore **-0.3
ignore 10.0-**
lmod relxill ~/data/monk/plot/warmcorona/model/relxill/relxill_model_v2.3
model zTBabs*TBabs(nthComp + atable{{{model_path}}} + xillver)
{param_lines}
{extra_cmds}
""".strip() + "\n"

    out_path.write_text(content, encoding="utf-8")
    return out_path


# =========================
# 8) 解析 show par（结构化，用于对齐 error）
# =========================
def parse_showpar_lines(lines: List[str]) -> List[dict]:
    """
    你要保留“原始 show par block”我们会另外原样打印；
    这里仅用于做一个结构化 recs，方便把 error(lo,hi) 对齐到所有参数上。
    """
    recs: List[dict] = []
    for line in lines:
        s = line.strip()
        if not s:
            continue

        # 常见 show par 的头部/分隔线/空行都跳过
        low = s.lower()
        if low.startswith("model ") or low.startswith("par") or low.startswith("model component") or low.startswith("source no"):
            continue
        if set(s) <= {"_", "-"}:
            continue

        toks = s.split()
        if len(toks) < 4:
            continue

        try:
            par_i = int(toks[0])
            comp_i = int(toks[1])
        except Exception:
            continue

        component = toks[2]

        # 从左到右找第一个 float 作为 value
        float_pos = None
        for j in range(3, len(toks)):
            if _FLOAT_RE.fullmatch(toks[j]):
                float_pos = j
                break
        if float_pos is None:
            continue

        parameter = " ".join(toks[3:float_pos]).strip()
        try:
            value = float(toks[float_pos])
        except Exception:
            continue

        tail_tokens = toks[float_pos + 1:]
        frozen = any(t.lower() == "frozen" for t in tail_tokens)

        link_to: Optional[str] = None
        # show par 通常是 "= p2"
        if "=" in tail_tokens:
            k = tail_tokens.index("=")
            if k + 1 < len(tail_tokens):
                link_to = tail_tokens[k + 1]

        pm: Optional[float] = None
        if "+/-" in tail_tokens:
            k = tail_tokens.index("+/-")
            if k + 1 < len(tail_tokens) and _FLOAT_RE.fullmatch(tail_tokens[k + 1]):
                try:
                    pm = float(tail_tokens[k + 1])
                except Exception:
                    pm = None

        recs.append({
            "par": par_i,
            "comp": comp_i,
            "component": component,
            "parameter": parameter,
            "value": value,
            "pm": pm,
            "frozen": frozen,
            "link": link_to,
            "raw": line.rstrip("\n"),
        })
    return recs


# =========================
# 9) 解析 XSPEC `error` 输出
# =========================
def parse_error_output(lines: List[str]) -> Dict[int, dict]:
    """
    尽量鲁棒地抓 error 输出。关键改进：
    - 先把 "= p2" / "p2" 这类 token 去掉，避免你之前看到的 err 里出现莫名其妙的 2/4
    - 只认“行首参数号 + 后面数字”
    """
    out: Dict[int, dict] = {}

    for line in lines:
        s0 = line.strip()
        if not s0:
            continue

        mpar = re.match(r"^\s*(\d+)\b", s0)
        if not mpar:
            continue
        par = int(mpar.group(1))

        # 清理掉 p2 之类
        s = re.sub(r"=\s*p\d+\b", "", s0, flags=re.I)
        s = re.sub(r"\bp\d+\b", "", s, flags=re.I)

        nums = _FLOAT_RE.findall(s)

        # 去掉第一个等于 par 的数
        if nums:
            try:
                if int(float(nums[0])) == par:
                    nums = nums[1:]
            except Exception:
                pass

        vals: List[float] = []
        for x in nums:
            try:
                vals.append(float(x))
            except Exception:
                vals = []
                break

        if len(vals) < 2:
            continue

        rec = {"kind": "unknown", "dminus": None, "dplus": None, "lo": None, "hi": None}

        # XSPEC error 常见两种：
        # (A) d- d+ lo hi
        # (B) lo hi
        if len(vals) >= 4:
            dminus, dplus, lo, hi = vals[0], vals[1], vals[2], vals[3]
            if lo > hi:
                # 兜底：有些输出顺序可能不同，按前两位当 lo/hi
                lo, hi, dminus, dplus = vals[0], vals[1], vals[2], vals[3]
                rec["kind"] = "4float(lo,hi,d-,d+)"
            else:
                rec["kind"] = "4float(d-,d+,lo,hi)"
            rec.update({"dminus": dminus, "dplus": dplus, "lo": lo, "hi": hi})
        else:
            lo, hi = vals[0], vals[1]
            rec.update({"lo": lo, "hi": hi, "kind": "2float(lo,hi)"})

        out[par] = rec

    return out


# =========================
# 10) XSPEC runner
# =========================
def _run_xspec_script(cwd: Path, script: str) -> Tuple[int, str]:
    proc = subprocess.run(
        ["xspec"],
        input=script,
        text=True,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return proc.returncode, proc.stdout


def _extract_block(raw: str, begin: str, end: str) -> List[str]:
    out: List[str] = []
    in_blk = False
    for line in raw.splitlines():
        if begin in line:
            in_blk = True
            continue
        if end in line:
            in_blk = False
            continue
        if in_blk:
            out.append(line.rstrip("\n"))
    return out


# =========================
# 11) 两次 XSPEC：Pass1 拿 free；Pass2 做 error + show par
# =========================
def run_xspec_pass1_free(
    source: str,
    obs_dir: Path,
    tom_xcm: Path,
    chi_xcm: Path,
    fit_iters: int,
) -> Tuple[List[int], int, str]:
    with tempfile.TemporaryDirectory(prefix="xmm_xspec_free_") as td:
        td = Path(td)

        tom_text = patch_text_for_mrk509(source, tom_xcm.read_text(encoding="utf-8", errors="ignore"))
        chi_text = patch_text_for_mrk509(source, chi_xcm.read_text(encoding="utf-8", errors="ignore"))

        tom2 = td / ("TOM_" + tom_xcm.name)
        chi2 = td / ("CHI_" + chi_xcm.name)
        tom2.write_text(tom_text, encoding="utf-8")
        chi2.write_text(chi_text, encoding="utf-8")

        tcl = f"""
query yes
@{tom2.as_posix()}
@{chi2.as_posix()}
fit {fit_iters}
puts "<<<PAR_TABLE_BEGIN>>>"
show par
puts "<<<PAR_TABLE_END>>>"
exit
""".lstrip()

        rc, raw = _run_xspec_script(obs_dir, tcl)

    par_lines = _extract_block(raw, "<<<PAR_TABLE_BEGIN>>>", "<<<PAR_TABLE_END>>>")
    recs = parse_showpar_lines(par_lines)
    free_list = [r["par"] for r in recs if (not r.get("frozen")) and (not r.get("link"))]
    free_list = sorted(set(free_list))
    return free_list, rc, raw


def run_xspec_pass2_error_and_showpar(
    source: str,
    obs_dir: Path,
    tom_xcm: Path,
    chi_xcm: Path,
    fit_iters: int,
    error_all: bool,
    free_list: List[int],
) -> Tuple[List[str], List[dict], Dict[int, dict], int, str, str]:
    """
    返回：
      - showpar_raw_block_lines（原样 show par）
      - showpar_recs（结构化）
      - errors（解析后的 error）
      - rc/raw
      - error_cmd_str
    """
    if error_all:
        error_cmd = "set n [tclout npar]\nerror 1-$n"
        error_cmd_str = "error 1-$npar (ALL)"
    else:
        if not free_list:
            error_cmd = "set n [tclout npar]\nerror 1-$n"
            error_cmd_str = "error 1-$npar (fallback, free empty)"
        else:
            parspec = compress_ranges(free_list)
            error_cmd = f"error {parspec}"
            error_cmd_str = f"error {parspec} (FREE)"

    with tempfile.TemporaryDirectory(prefix="xmm_xspec_err_") as td:
        td = Path(td)

        tom_text = patch_text_for_mrk509(source, tom_xcm.read_text(encoding="utf-8", errors="ignore"))
        chi_text = patch_text_for_mrk509(source, tom_xcm.read_text(encoding="utf-8", errors="ignore"))  # placeholder
        chi_text = patch_text_for_mrk509(source, chi_xcm.read_text(encoding="utf-8", errors="ignore"))

        tom2 = td / ("TOM_" + tom_xcm.name)
        chi2 = td / ("CHI_" + chi_xcm.name)
        tom2.write_text(tom_text, encoding="utf-8")
        chi2.write_text(chi_text, encoding="utf-8")

        tcl = f"""
query yes
@{tom2.as_posix()}
@{chi2.as_posix()}
fit {fit_iters}
puts "<<<ERROR_OUTPUT_BEGIN>>>"
{error_cmd}
puts "<<<ERROR_OUTPUT_END>>>"
puts "<<<SHOWPAR_RAW_BEGIN>>>"
show par
puts "<<<SHOWPAR_RAW_END>>>"
exit
""".lstrip()

        rc, raw = _run_xspec_script(obs_dir, tcl)

    # 原样 block
    showpar_raw = _extract_block(raw, "<<<SHOWPAR_RAW_BEGIN>>>", "<<<SHOWPAR_RAW_END>>>")

    # 结构化 recs：从原样 block 里解析（用于对齐 error）
    recs = parse_showpar_lines(showpar_raw)

    # error block
    err_lines = _extract_block(raw, "<<<ERROR_OUTPUT_BEGIN>>>", "<<<ERROR_OUTPUT_END>>>")
    errors = parse_error_output(err_lines)

    return showpar_raw, recs, errors, rc, raw, error_cmd_str


# =========================
# 12) 目标解析：source obsid obsid source obsid...
# =========================
def parse_targets(tokens: List[str]) -> Dict[str, List[str]]:
    targets: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    for t in tokens:
        if _OBSID_RE.match(t):
            if cur is None:
                raise ValueError(f"Found OBSID '{t}' before any SOURCE.")
            targets.setdefault(cur, []).append(t)
        else:
            cur = t
            targets.setdefault(cur, [])
    return targets


# =========================
# 13) 输出：对齐所有参数的 error
# =========================
def format_all_params_with_error(recs: List[dict], errors: Dict[int, dict]) -> List[str]:
    lines = []
    lines.append(f"{'par':>3s} {'comp':>4s}  {'component':<18s} {'parameter':<24s} {'value':>14s}  {'pm':>14s}  {'err(lo,hi)':>22s}  flags")
    for r in sorted(recs, key=lambda x: x["par"]):
        par_i = r["par"]
        val = r["value"]
        pm = r["pm"]
        err = errors.get(par_i)

        sval = f"{val:.6g}"
        spm = "--" if pm is None else f"{pm:.6g}"
        serr = "--"
        if err is not None and err.get("lo") is not None and err.get("hi") is not None:
            serr = f"{err['lo']:.6g},{err['hi']:.6g}"

        flags = []
        if r.get("frozen"):
            flags.append("frozen")
        if r.get("link"):
            flags.append(f"= {r['link']}")
        if err is not None:
            flags.append(f"err:{err.get('kind','?')}")

        lines.append(
            f"{par_i:3d} {r['comp']:4d}  {r['component'][:18]:<18s} {r['parameter'][:24]:<24s} "
            f"{sval:>14s}  {spm:>14s}  {serr:>22s}  {' '.join(flags)}"
        )
    return lines


# =========================
# 14) main
# =========================
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Scan XMM dir, pick chi_*.xcm (prefer 1 sphere+1 slab; exclude gdiskwien), "
            "generate Tom xcm via online z/NH, run XSPEC.\n"
            "If --do-error: TWO-PASS XSPEC.\n"
            "Output keeps RAW show par block + then PASS summary + then all-parameter errors."
        )
    )

    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    ap.add_argument("--recursive", action="store_true", help="递归搜 chi_*.xcm（默认只搜 obsid/ 和 obsid/xcm/）")

    ap.add_argument("--read-xcm", action="store_true", help="读 chi_*.xcm 并跑 XSPEC dump（否则只做 index）")
    ap.add_argument("--all", action="store_true", help="对所有 source/obsid 执行（慎用）")

    ap.add_argument("--debug-xspec-out", action="store_true", help="保存完整 XSPEC log（排错用）")
    ap.add_argument("--include-cpd", action="store_true", help="写 Tom xcm 时包含 'cpd /xs'")

    ap.add_argument("--do-fit", action="store_true", help="在 @chi 后执行 fit（若 --do-error 则会强制 fit）")
    ap.add_argument("--do-error", action="store_true", help="执行 error（默认：ALL params 1-$npar）")
    ap.add_argument("--error-free-only", action="store_true", help="改为只对 free & unlinked 的参数跑 error（你以前的模式）")
    ap.add_argument("--fit-iters", type=int, default=1000, help="fit 最大迭代次数（--do-error 会强制 fit）")

    ap.add_argument("--no-online", action="store_true", help="禁用在线查 z/NH（将使用 z=0, nh=0.01）")

    ap.add_argument("tokens", nargs="*", help="如：1H0419-577 0148000201 0148000401  Mrk509 0130720101")

    args = ap.parse_args()

    out_root = args.out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    # (A) 扫描 + 写 index（始终生成）
    entries = scan_all_entries()
    full_idx = build_full_index(entries, recursive=args.recursive)
    index_txt = write_full_index_txt(out_root, entries, full_idx)
    print(f"[OK] wrote full index: {index_txt}")

    if not args.read_xcm:
        return

    # (B) 选 targets
    if args.all:
        selected = entries
    else:
        targets = parse_targets(args.tokens)
        by_source: Dict[str, List[ObsEntry]] = {}
        for e in entries:
            by_source.setdefault(e.source, []).append(e)

        selected = []
        for src_in, obsids in targets.items():
            real_src = None
            for k in by_source.keys():
                if k.lower() == src_in.lower():
                    real_src = k
                    break
            if real_src is None:
                continue
            if obsids:
                sset = set(obsids)
                selected.extend([e for e in by_source[real_src] if e.obsid in sset])
            else:
                selected.extend(by_source[real_src])

    uniq = {(e.source, e.obsid): e for e in selected}
    selected = [uniq[k] for k in sorted(uniq.keys())]

    dump_base = out_root / "xmm_paramdump_with_error"
    dump_base.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    srcs = sorted(set(e.source for e in selected))
    src_tag = srcs[0] if len(srcs) == 1 else f"{len(srcs)}src"
    obs_tag = f"{len(selected)}obs"
    dump_txt = dump_base / f"xspec_paramdump_with_error_{src_tag}_{obs_tag}_{stamp}.txt"

    debug_dir = dump_base / "debug_xspec"
    if args.debug_xspec_out:
        debug_dir.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append("XSPEC param dump: keep RAW show par + PASS summary + ALL-param error table")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Selected targets: {len(selected)}")
    do_fit_effective = args.do_fit or args.do_error
    lines.append(f"do_fit={do_fit_effective}  do_error={args.do_error}  fit_iters={args.fit_iters}")
    lines.append(f"online_znh={'OFF' if args.no_online else 'ON'}")
    if args.do_error:
        lines.append(f"error_mode: {'FREE-only' if args.error_free_only else 'ALL(1-$npar)'}")
    lines.append("")

    # 每源缓存 z/nh
    znh_cache: Dict[str, Tuple[float, float, str, str]] = {}
    total_tables = 0

    for e in selected:
        picks, note = full_idx.get((e.source, e.obsid), ([], ""))
        if not picks:
            continue

        spec = e.path / "PN_spectrum_grp.fits"
        if not spec.exists():
            continue

        # z/nh cache
        if e.source not in znh_cache:
            if args.no_online:
                znh_cache[e.source] = (0.0, 0.01, e.source, "")
            else:
                try:
                    sim = query_simbad_basic(e.source)
                    z = float(sim["z"]) if sim["z"] is not None else 0.0
                    nh = fetch_w3nh_nh(float(sim["ra_deg"]), float(sim["dec_deg"]), radius_deg=0.1)
                    nh_used = nh["wavg_nh"] if nh["wavg_nh"] is not None else nh["avg_nh"]
                    if nh_used is None:
                        nh_used = 1.0e20
                    nh_1e22 = float(nh_used) / 1e22
                    znh_cache[e.source] = (z, nh_1e22, sim.get("used_name", e.source), nh["url"])
                except Exception as ex:
                    znh_cache[e.source] = (0.0, 0.01, e.source, "")
                    lines.append(f"[WARN] z/NH query failed for {e.source}: {ex} -> fallback z=0, nh=0.01")

        z, nh_1e22, sim_used_name, w3nh_url = znh_cache[e.source]
        params = build_params_for_source(z, nh_1e22)

        # 写 Tom xcm（slab/sphere）
        tom_slab = write_tom_xcm(e.path, e.source, e.obsid, slab=True, params=params, include_cpd=args.include_cpd)
        tom_sph = write_tom_xcm(e.path, e.source, e.obsid, slab=False, params=params, include_cpd=args.include_cpd)

        for rank, pick in enumerate(sorted(picks, key=lambda x: x.mtime, reverse=True), 1):
            kind = pick.kind
            tom = tom_sph if kind == "sphere" else tom_slab

            free_list: List[int] = []
            rc1, raw1 = 0, ""

            if args.do_error:
                # Pass1：只为 FREE-only 模式准备 free_list；ALL 模式也保留 free_list 方便你记录
                free_list, rc1, raw1 = run_xspec_pass1_free(
                    source=e.source,
                    obs_dir=e.path,
                    tom_xcm=tom,
                    chi_xcm=pick.path,
                    fit_iters=args.fit_iters,
                )

                # Pass2：error + raw show par
                showpar_raw, recs, errors, rc2, raw2, errcmd = run_xspec_pass2_error_and_showpar(
                    source=e.source,
                    obs_dir=e.path,
                    tom_xcm=tom,
                    chi_xcm=pick.path,
                    fit_iters=args.fit_iters,
                    error_all=(not args.error_free_only),
                    free_list=free_list,
                )
                rc = rc2
                raw = raw2
            else:
                # 不做 error：仍然想要你要的 show par 原样输出
                # 用 Pass2 的框架，但把 error_cmd 变成空：这里简单起见，直接执行一次：fit? + show par
                # （保持结构最少改动）
                with tempfile.TemporaryDirectory(prefix="xmm_xspec_nerr_") as td:
                    td = Path(td)
                    tom_text = patch_text_for_mrk509(e.source, tom.read_text(encoding="utf-8", errors="ignore"))
                    chi_text = patch_text_for_mrk509(e.source, pick.path.read_text(encoding="utf-8", errors="ignore"))
                    tom2 = td / ("TOM_" + tom.name)
                    chi2 = td / ("CHI_" + pick.path.name)
                    tom2.write_text(tom_text, encoding="utf-8")
                    chi2.write_text(chi_text, encoding="utf-8")
                    fit_line = f"fit {args.fit_iters}" if args.do_fit else ""
                    tcl = f"""
query yes
@{tom2.as_posix()}
@{chi2.as_posix()}
{fit_line}
puts "<<<SHOWPAR_RAW_BEGIN>>>"
show par
puts "<<<SHOWPAR_RAW_END>>>"
exit
""".lstrip()
                    rc, raw = _run_xspec_script(e.path, tcl)
                showpar_raw = _extract_block(raw, "<<<SHOWPAR_RAW_BEGIN>>>", "<<<SHOWPAR_RAW_END>>>")
                recs = parse_showpar_lines(showpar_raw)
                errors = {}
                errcmd = "(none)"

            if args.debug_xspec_out:
                (debug_dir / f"xspec_{e.source}_{e.obsid}_rank{rank}_{pick.path.name}.log").write_text(raw, encoding="utf-8")
                if args.do_error:
                    (debug_dir / f"xspec_{e.source}_{e.obsid}_rank{rank}_{pick.path.name}.pass1_free.log").write_text(raw1, encoding="utf-8")

            lines.append("=" * 112)
            lines.append(f"{e.source}/{e.obsid}  rank={rank}  chi_kind={kind}  chi_mtime={pick.mtime_str}")
            lines.append(f"pick_note: {note}")
            lines.append(f"TOM: {tom}")
            lines.append(f"CHI: {pick.path}")
            lines.append(f"SIMBAD_used_name: {sim_used_name}   z={z:.6g}")
            lines.append(f"w3nh_url: {w3nh_url}   nh(TBabs)={nh_1e22:.6g} (1e22)")
            lines.append("-" * 112)

            # ① 你要保留的：RAW show par（原样）
            lines.append("<<<SHOWPAR_RAW_BEGIN>>>")
            if showpar_raw:
                lines.extend(showpar_raw)
            else:
                lines.append("(EMPTY SHOWPAR BLOCK)")
            lines.append("<<<SHOWPAR_RAW_END>>>")
            lines.append("")

            # ② 再输出你要的 PASS 信息 + short 摘要
            if args.do_error:
                lines.append(f"PASS1(xspec) rc={rc1}  free_count={len(free_list)}  free_pars={compress_ranges(free_list) if free_list else '(EMPTY)'}")
                lines.append(f"PASS2(xspec) rc={rc}   error_cmd={errcmd}")
            else:
                lines.append(f"PASS(xspec) rc={rc}  error_cmd={errcmd}")

            fitinfo = parse_fit_summary(raw)
            lines.append(f"fit: chi2={fitinfo['chi2']}  dof={fitinfo['dof']}  nhp={fitinfo['nhp']}  bins={fitinfo['bins']}")
            lines.append(f"free_pars: {compress_ranges(free_list) if free_list else '(n/a)'}")
            lines.append("")

            # ③ 你现在要的：所有参数的 error（对齐到 show par recs）
            if recs:
                total_tables += 1
                lines.append("params(all): show par table + error(lo,hi) for ALL parameters (if available)")
                lines.extend(format_all_params_with_error(recs, errors))
            else:
                lines.append("(NO structured params parsed from show par. But RAW show par is kept above.)")

            lines.append("")

    dump_txt.write_text("\n".join(lines), encoding="utf-8")
    print(f"[OK] wrote param dump: {dump_txt}")
    print(f"[INFO] parsed show-par tables: {total_tables}")


if __name__ == "__main__":
    main()
