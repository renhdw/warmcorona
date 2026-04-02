#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xmm_index_v2.py

功能概述
--------
这是一个用于 XMM-Newton 光谱拟合批处理的自动化脚本，主要完成以下工作：

1. 扫描 XMM 数据目录
2. 查找 chi_*.xcm 文件（优先 sphere + slab）
3. 自动查询 redshift / Galactic NH（可关闭）
4. 自动生成 Tom 模型 xcm
5. 自动调用 XSPEC 执行 fit / error
6. 输出索引、txt 参数汇总、观测级 summary csv
7. 支持并行 (--jobs)
8. 支持进度追踪 (--track-progress)

进度追踪功能
------------
开启 --track-progress 后：
- 终端实时显示任务开始 / pass1 / pass2 / 完成
- 在 debug_xspec/ 下写入每个任务的状态文件
- 写入总进度汇总文件 progress_summary.txt

常用命令
--------
1) 只扫描:
   python3 xmm_index_v2.py

2) 单目标:
   python3 xmm_index_v2.py --read-xcm --do-fit --do-error --error-free-only Mrk509 0130720101

3) 并行:
   python3 xmm_index_v2.py --read-xcm --do-fit --do-error --error-free-only --jobs 4 Mrk509 0130720101

4) 跑全部并追踪进度:
   python3 xmm_index_v2.py --read-xcm --do-fit --do-error --error-free-only --all --jobs 6 --track-progress

推荐
----
你的机器（20 logical CPU）建议：
- 单个命令: --jobs 4 或 6
- 多个命令同时跑: 每个命令 --jobs 1 或 2
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import argparse
import re
import subprocess
import tempfile
from datetime import datetime
import urllib.parse as up
import csv
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import requests
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u


# =========================
# 0) 环境路径（按需改）
# =========================
XMM_BASE_DIR = Path("/home/hdw/data/XMM").resolve()
DEFAULT_OUT_ROOT = Path("~/data/monk/plot/warmcorona/test/2026.04.01").expanduser().resolve()

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

_EXCLUDE_CHI_SUBSTRS = ["gdiskwien"]

_FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
_SCI_RE = re.compile(r"([0-9.+-]+E[+-]?[0-9]+)", re.IGNORECASE)

_FITSTAT_RE = re.compile(r"Fit statistic\s*:\s*Chi-Squared\s+([0-9.]+)\s+using\s+(\d+)\s+bins", re.I)
_DOF_RE = re.compile(r"Null hypothesis probability of\s+([0-9.eE+-]+)\s+with\s+(\d+)\s+degrees of freedom", re.I)

MODEL_LINE_RE = re.compile(r"^Model\s+(.+?)\s+Source No\.:", re.I)


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


def safe_float_str(x: Any) -> str:
    if x is None:
        return ""
    try:
        return f"{float(x):.12g}"
    except Exception:
        return str(x)


def bool01(x: bool) -> int:
    return 1 if x else 0


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# =========================
# 4) 进度追踪 helpers
# =========================
def task_tag(source: str, obsid: str, geom: str) -> str:
    return f"{source}_{obsid}_{geom}"


def write_status_file(debug_dir: Optional[Path], source: str, obsid: str, geom: str, text: str) -> None:
    if debug_dir is None:
        return
    debug_dir.mkdir(parents=True, exist_ok=True)
    p = debug_dir / f"status_{task_tag(source, obsid, geom)}.txt"
    p.write_text(text, encoding="utf-8")


def append_status_file(debug_dir: Optional[Path], source: str, obsid: str, geom: str, text: str) -> None:
    if debug_dir is None:
        return
    debug_dir.mkdir(parents=True, exist_ok=True)
    p = debug_dir / f"status_{task_tag(source, obsid, geom)}.txt"
    with p.open("a", encoding="utf-8") as f:
        f.write(text)


def write_progress_summary(
    debug_dir: Optional[Path],
    total_tasks: int,
    done_tasks: int,
    running_tasks: int,
    failed_tasks: int,
    finished_tags: List[str],
) -> None:
    if debug_dir is None:
        return
    debug_dir.mkdir(parents=True, exist_ok=True)
    p = debug_dir / "progress_summary.txt"
    lines = [
        f"time = {now_str()}",
        f"total_tasks = {total_tasks}",
        f"done_tasks = {done_tasks}",
        f"running_tasks = {running_tasks}",
        f"failed_tasks = {failed_tasks}",
        "",
        "finished:",
    ]
    lines.extend(finished_tags[-50:])
    p.write_text("\n".join(lines), encoding="utf-8")


# =========================
# 5) 扫描
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
# 6) 挑 chi_xcm
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
# 7) 在线查 z/NH
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
# 8) 写 Tom xcm
# =========================
def build_params_for_source(z: Optional[float], nh_1e22: Optional[float]) -> List[float]:
    if z is None:
        z = 0.0
    if nh_1e22 is None:
        nh_1e22 = 0.01
    return [
        1.0, z, nh_1e22,
        1.5, 100.0, 3e-3, 0,
        z, 1.0,
        0.5, 12.0, z, 1.0,
        1.5, 1.0, 300.0, 0.0, z, 30.0, -1.0, 1.0
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
# 9) 解析 show par
# =========================
def parse_showpar_lines(lines: List[str]) -> List[dict]:
    recs: List[dict] = []
    for line in lines:
        s = line.strip()
        if not s:
            continue

        low = s.lower()
        if low.startswith("model ") or low.startswith("par") or low.startswith("model component") or low.startswith("source no"):
            continue
        if set(s) <= {"_", "-", "="}:
            continue
        if low.startswith("xspec12>show par"):
            continue
        if low.startswith("parameters defined"):
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


def parse_model_line(showpar_raw: List[str]) -> str:
    for line in showpar_raw:
        m = MODEL_LINE_RE.search(line.strip())
        if m:
            return m.group(1).strip()
    return ""


# =========================
# 10) 解析 XSPEC error
# =========================
def parse_error_output(lines: List[str]) -> Dict[int, dict]:
    out: Dict[int, dict] = {}

    for line in lines:
        s0 = line.strip()
        if not s0:
            continue

        mpar = re.match(r"^\s*(\d+)\b", s0)
        if not mpar:
            continue
        par = int(mpar.group(1))

        s = re.sub(r"=\s*p\d+\b", "", s0, flags=re.I)
        s = re.sub(r"\bp\d+\b", "", s, flags=re.I)

        nums = _FLOAT_RE.findall(s)
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

        if len(vals) >= 4:
            dminus, dplus, lo, hi = vals[0], vals[1], vals[2], vals[3]
            if lo > hi:
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


def err_to_str(err: Optional[dict]) -> str:
    if not err:
        return ""
    lo = err.get("lo")
    hi = err.get("hi")
    if lo is None or hi is None:
        return ""
    return f"{safe_float_str(lo)},{safe_float_str(hi)}"


# =========================
# 11) XSPEC runner
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
# 12) 两次 XSPEC
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

    showpar_raw = _extract_block(raw, "<<<SHOWPAR_RAW_BEGIN>>>", "<<<SHOWPAR_RAW_END>>>")
    recs = parse_showpar_lines(showpar_raw)
    err_lines = _extract_block(raw, "<<<ERROR_OUTPUT_BEGIN>>>", "<<<ERROR_OUTPUT_END>>>")
    errors = parse_error_output(err_lines)

    return showpar_raw, recs, errors, rc, raw, error_cmd_str


# =========================
# 13) 目标解析
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
# 14) txt 输出
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
# 15) summary helpers
# =========================
def normalize_component_name(component: str) -> str:
    c = (component or "").strip()
    if c.startswith("warmcom_raw"):
        return "warmcom"
    if c.startswith("warmcom_pure"):
        return "warmcom"
    return c


def warmcom_mode_from_component(component: str) -> str:
    c = (component or "").strip()
    if c.startswith("warmcom_raw"):
        return "raw"
    if c.startswith("warmcom_pure"):
        return "pure"
    return ""


def simplify_param_key(component: str, parameter: str) -> Optional[str]:
    comp = normalize_component_name(component)
    p = (parameter or "").strip().lower()

    if comp == "zTBabs":
        if "nh" in p:
            return "zt_nh"
        if "redshift" in p:
            return "zt_z"

    if comp == "TBabs":
        if "nh" in p:
            return "tb_nh"

    if comp == "nthComp":
        if p == "gamma":
            return "nth_gam"
        if "kt_e" in p or "kt e" in p:
            return "nth_kte"
        if "kt_bb" in p or "kt bb" in p:
            return "nth_kbb"
        if "inp_type" in p:
            return "nth_inp"
        if p == "redshift":
            return "nth_z"
        if p == "norm":
            return "nth_norm"

    if comp == "warmcom":
        if p == "":
            if component.endswith("te"):
                return "wc_te"
            if component.endswith("tau"):
                return "wc_tau"
            if component.endswith("z"):
                return "wc_z"
            if component.endswith("norm"):
                return "wc_norm"
        if "te" in p:
            return "wc_te"
        if "tau" in p:
            return "wc_tau"
        if p == "z" or "redshift" in p:
            return "wc_z"
        if "norm" in p:
            return "wc_norm"

    if comp == "zgauss":
        if "linee" in p:
            return "g_e"
        if "sigma" in p:
            return "g_sig"
        if "redshift" in p:
            return "g_z"
        if p == "norm":
            return "g_norm"

    if comp == "xillver":
        if p == "gamma":
            return "x_gam"
        if p == "afe":
            return "x_afe"
        if "ecut" in p:
            return "x_ecut"
        if "logxi" in p:
            return "x_xi"
        if p == "z":
            return "x_z"
        if "incl" in p:
            return "x_inc"
        if "refl_frac" in p:
            return "x_refl"
        if p == "norm":
            return "x_norm"

    return None


def model_flags_from_recs(recs: List[dict], model_expr: str, chi_path: Path) -> Dict[str, Any]:
    comps = [normalize_component_name(r.get("component", "")) for r in recs]
    comps_lower = [c.lower() for c in comps]
    expr = (model_expr or "").lower()
    chi_text = str(chi_path).lower()

    has_gauss = ("zgauss" in comps_lower) or ("zgauss" in expr)
    has_xillver = ("xillver" in comps_lower) or ("xillver" in expr)
    has_nth = ("nthcomp" in [c.lower() for c in comps]) or ("nthcomp" in expr)
    has_warm = any(c == "warmcom" for c in comps) or ("atable{" in expr)
    has_wa = ("wa_" in chi_text) or ("wa_" in expr) or ("cloudy" in chi_text)

    wa_name = ""
    chi_name = chi_path.name
    if "wa_Mrk509_xi_NH.mod".lower() in chi_name.lower():
        wa_name = "wa_Mrk509_xi_NH"
    elif "wa_cloudy_agn_v2wa_agn_xi_NH.mod".lower() in chi_name.lower():
        wa_name = "wa_cloudy_agn_v2wa_agn_xi_NH"
    elif "wa_" in chi_name.lower():
        wa_name = chi_name

    return {
        "has_wa": bool01(has_wa),
        "has_gauss": bool01(has_gauss),
        "has_xillver": bool01(has_xillver),
        "has_nth": bool01(has_nth),
        "has_warm": bool01(has_warm),
        "wa_name": wa_name,
    }


def build_model_label(geom: str, flags: Dict[str, Any]) -> str:
    parts = [geom if geom else "unknown"]
    if flags.get("has_wa"):
        parts.append("wa")
    if flags.get("has_nth"):
        parts.append("nth")
    if flags.get("has_warm"):
        parts.append("warm")
    if flags.get("has_gauss"):
        parts.append("gauss")
    if flags.get("has_xillver"):
        parts.append("xill")
    return "+".join(parts)


def ensure_param_columns(row: Dict[str, Any], base: str) -> None:
    row.setdefault(base, "")
    row.setdefault(f"{base}_free", "")
    row.setdefault(f"{base}_err", "")
    row.setdefault(f"{base}_pm", "")


def fill_param_row(row: Dict[str, Any], base: str, rec: dict, err: Optional[dict]) -> None:
    ensure_param_columns(row, base)
    row[base] = rec.get("value", "")
    row[f"{base}_free"] = 0 if (rec.get("frozen") or rec.get("link")) else 1
    row[f"{base}_err"] = err_to_str(err)
    row[f"{base}_pm"] = rec.get("pm", "")


def build_empty_obs_row(
    source: str,
    obsid: str,
    geom: str,
    z: float,
    nh_1e22: float,
    simbad_used_name: str,
    pick_note: str,
) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        "source": source,
        "obsid": obsid,
        "geom": geom,
        "rank": "",
        "chi_kind": geom,
        "chi": "",
        "dof": "",
        "chi_dof": "",
        "nhp": "",
        "bins": "",
        "z": z,
        "nh_gal": nh_1e22,
        "model": "",
        "model_expr": "",
        "simbad_name": simbad_used_name,
        "pick_note": pick_note,
        "chi_path": "",
        "tom_path": "",
        "chi_mtime": "",
        "free_pars": "",
        "free_n": "",
        "has_wa": "",
        "has_gauss": "",
        "has_xillver": "",
        "has_nth": "",
        "has_warm": "",
        "wa_name": "",
        "wc_mode": "",
    }

    short_cols = [
        "zt_nh", "zt_z", "tb_nh",
        "nth_gam", "nth_kte", "nth_kbb", "nth_inp", "nth_z", "nth_norm",
        "wc_te", "wc_tau", "wc_z", "wc_norm",
        "g_e", "g_sig", "g_z", "g_norm",
        "x_gam", "x_afe", "x_ecut", "x_xi", "x_z", "x_inc", "x_refl", "x_norm",
    ]
    for b in short_cols:
        ensure_param_columns(row, b)

    return row


def build_obs_summary_row(
    source: str,
    obsid: str,
    geom: str,
    rank: int,
    chi_kind: str,
    chi_path: Path,
    tom_path: Path,
    chi_mtime_str: str,
    pick_note: str,
    simbad_used_name: str,
    z: float,
    nh_1e22: float,
    fitinfo: Dict[str, object],
    free_list: List[int],
    recs: List[dict],
    errors: Dict[int, dict],
    model_expr: str,
) -> Dict[str, Any]:
    row = build_empty_obs_row(source, obsid, geom, z, nh_1e22, simbad_used_name, pick_note)

    row["rank"] = rank
    row["chi_kind"] = chi_kind
    row["chi"] = fitinfo.get("chi2", "")
    row["dof"] = fitinfo.get("dof", "")
    row["nhp"] = fitinfo.get("nhp", "")
    row["bins"] = fitinfo.get("bins", "")
    if fitinfo.get("chi2") is not None and fitinfo.get("dof") not in (None, 0):
        row["chi_dof"] = float(fitinfo["chi2"]) / float(fitinfo["dof"])
    else:
        row["chi_dof"] = ""

    row["model_expr"] = model_expr
    row["simbad_name"] = simbad_used_name
    row["chi_path"] = str(chi_path)
    row["tom_path"] = str(tom_path)
    row["chi_mtime"] = chi_mtime_str
    row["free_pars"] = compress_ranges(free_list) if free_list else ""
    row["free_n"] = len(free_list)

    flags = model_flags_from_recs(recs, model_expr, chi_path)
    row.update(flags)
    row["model"] = build_model_label(geom, flags)

    for r in recs:
        base = simplify_param_key(r.get("component", ""), r.get("parameter", ""))
        if base is None:
            base = simplify_param_key(r.get("component", ""), "")
        if base is None:
            continue
        fill_param_row(row, base, r, errors.get(r["par"]))

        wc_mode = warmcom_mode_from_component(r.get("component", ""))
        if wc_mode:
            row["wc_mode"] = wc_mode

    return row


def write_csv_rows(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames: List[str] = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)

    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            clean = {}
            for k in fieldnames:
                v = row.get(k, "")
                if isinstance(v, float):
                    clean[k] = safe_float_str(v)
                elif v is None:
                    clean[k] = ""
                else:
                    clean[k] = v
            w.writerow(clean)


# =========================
# 16) 单任务 worker
# =========================
def worker_run_fit_task(task: Dict[str, Any]) -> Dict[str, Any]:
    source = task["source"]
    obsid = task["obsid"]
    geom = task["geom"]
    rank = task["rank"]
    chi_kind = task["chi_kind"]
    obs_dir = Path(task["obs_dir"])
    chi_path = Path(task["chi_path"])
    tom_path = Path(task["tom_path"])
    note = task["pick_note"]
    sim_used_name = task["simbad_used_name"]
    z = task["z"]
    nh_1e22 = task["nh_1e22"]
    w3nh_url = task["w3nh_url"]
    do_error = task["do_error"]
    do_fit = task["do_fit"]
    error_free_only = task["error_free_only"]
    fit_iters = task["fit_iters"]
    debug_xspec_out = task["debug_xspec_out"]
    track_progress = task["track_progress"]
    debug_dir = Path(task["debug_dir"]) if task["debug_dir"] else None

    def log_state(msg: str) -> None:
        if track_progress:
            append_status_file(debug_dir, source, obsid, geom, f"[{now_str()}] {msg}\n")

    if track_progress:
        write_status_file(
            debug_dir, source, obsid, geom,
            f"[{now_str()}] START {source} {obsid} {geom}\n"
        )

    free_list: List[int] = []
    rc1, raw1 = 0, ""

    try:
        if do_error:
            log_state("PASS1_BEGIN fit + show par")
            free_list, rc1, raw1 = run_xspec_pass1_free(
                source=source,
                obs_dir=obs_dir,
                tom_xcm=tom_path,
                chi_xcm=chi_path,
                fit_iters=fit_iters,
            )
            log_state("PASS1_END")

            log_state("PASS2_BEGIN fit + error + show par")
            showpar_raw, recs, errors, rc2, raw2, errcmd = run_xspec_pass2_error_and_showpar(
                source=source,
                obs_dir=obs_dir,
                tom_xcm=tom_path,
                chi_xcm=chi_path,
                fit_iters=fit_iters,
                error_all=(not error_free_only),
                free_list=free_list,
            )
            log_state("PASS2_END")
            rc = rc2
            raw = raw2
        else:
            log_state("PASS_BEGIN fit/show par")
            with tempfile.TemporaryDirectory(prefix="xmm_xspec_nerr_") as td:
                td = Path(td)
                tom_text = patch_text_for_mrk509(source, tom_path.read_text(encoding="utf-8", errors="ignore"))
                chi_text = patch_text_for_mrk509(source, chi_path.read_text(encoding="utf-8", errors="ignore"))
                tom2 = td / ("TOM_" + tom_path.name)
                chi2 = td / ("CHI_" + chi_path.name)
                tom2.write_text(tom_text, encoding="utf-8")
                chi2.write_text(chi_text, encoding="utf-8")
                fit_line = f"fit {fit_iters}" if do_fit else ""
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
                rc, raw = _run_xspec_script(obs_dir, tcl)
            showpar_raw = _extract_block(raw, "<<<SHOWPAR_RAW_BEGIN>>>", "<<<SHOWPAR_RAW_END>>>")
            recs = parse_showpar_lines(showpar_raw)
            errors = {}
            errcmd = "(none)"
            log_state("PASS_END")

        if debug_xspec_out and debug_dir is not None:
            debug_dir.mkdir(parents=True, exist_ok=True)
            (debug_dir / f"xspec_{source}_{obsid}_{geom}_{chi_path.name}.log").write_text(raw, encoding="utf-8")
            if do_error:
                (debug_dir / f"xspec_{source}_{obsid}_{geom}_{chi_path.name}.pass1_free.log").write_text(raw1, encoding="utf-8")

        fitinfo = parse_fit_summary(raw)
        model_expr = parse_model_line(showpar_raw)

        txt_lines: List[str] = []
        txt_lines.append("=" * 112)
        txt_lines.append(f"{source}/{obsid}  rank={rank}  chi_kind={chi_kind}  chi_mtime={task['chi_mtime_str']}")
        txt_lines.append(f"pick_note: {note}")
        txt_lines.append(f"TOM: {tom_path}")
        txt_lines.append(f"CHI: {chi_path}")
        txt_lines.append(f"SIMBAD_used_name: {sim_used_name}   z={z:.6g}")
        txt_lines.append(f"w3nh_url: {w3nh_url}   nh(TBabs)={nh_1e22:.6g} (1e22)")
        txt_lines.append("-" * 112)
        txt_lines.append("<<<SHOWPAR_RAW_BEGIN>>>")
        if showpar_raw:
            txt_lines.extend(showpar_raw)
        else:
            txt_lines.append("(EMPTY SHOWPAR BLOCK)")
        txt_lines.append("<<<SHOWPAR_RAW_END>>>")
        txt_lines.append("")

        if do_error:
            txt_lines.append(f"PASS1(xspec) rc={rc1}  free_count={len(free_list)}  free_pars={compress_ranges(free_list) if free_list else '(EMPTY)'}")
            txt_lines.append(f"PASS2(xspec) rc={rc}   error_cmd={errcmd}")
        else:
            txt_lines.append(f"PASS(xspec) rc={rc}  error_cmd={errcmd}")

        txt_lines.append(f"fit: chi2={fitinfo['chi2']}  dof={fitinfo['dof']}  nhp={fitinfo['nhp']}  bins={fitinfo['bins']}")
        txt_lines.append(f"free_pars: {compress_ranges(free_list) if free_list else '(n/a)'}")
        txt_lines.append("")

        if recs:
            txt_lines.append("params(all): show par table + error(lo,hi) for ALL parameters (if available)")
            txt_lines.extend(format_all_params_with_error(recs, errors))
        else:
            txt_lines.append("(NO structured params parsed from show par. But RAW show par is kept above.)")
        txt_lines.append("")

        obs_row = build_obs_summary_row(
            source=source,
            obsid=obsid,
            geom=geom,
            rank=rank,
            chi_kind=chi_kind,
            chi_path=chi_path,
            tom_path=tom_path,
            chi_mtime_str=task["chi_mtime_str"],
            pick_note=note,
            simbad_used_name=sim_used_name,
            z=z,
            nh_1e22=nh_1e22,
            fitinfo=fitinfo,
            free_list=free_list,
            recs=recs,
            errors=errors,
            model_expr=model_expr,
        )

        log_state(f"DONE chi2={fitinfo.get('chi2')} dof={fitinfo.get('dof')}")

        return {
            "ok": True,
            "order": task["order"],
            "source": source,
            "obsid": obsid,
            "geom": geom,
            "txt_lines": txt_lines,
            "obs_row": obs_row,
            "has_table": bool(recs),
        }

    except Exception as ex:
        log_state(f"ERROR {ex}")
        err_lines = [
            "=" * 112,
            f"{source}/{obsid}  rank={rank}  chi_kind={chi_kind}",
            f"[ERROR] task failed: {ex}",
            ""
        ]
        return {
            "ok": False,
            "order": task["order"],
            "source": source,
            "obsid": obsid,
            "geom": geom,
            "txt_lines": err_lines,
            "obs_row": build_empty_obs_row(
                source=source,
                obsid=obsid,
                geom=geom,
                z=z,
                nh_1e22=nh_1e22,
                simbad_used_name=sim_used_name,
                pick_note=f"{note}; ERROR={ex}",
            ),
            "has_table": False,
        }


# =========================
# 17) main
# =========================
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Scan XMM dir, pick chi_*.xcm, generate Tom xcm via online z/NH, "
            "run XSPEC, support parallel jobs and optional progress tracking."
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
    ap.add_argument("--error-free-only", action="store_true", help="改为只对 free & unlinked 的参数跑 error")
    ap.add_argument("--fit-iters", type=int, default=1000, help="fit 最大迭代次数（--do-error 会强制 fit）")

    ap.add_argument("--no-online", action="store_true", help="禁用在线查 z/NH（将使用 z=0, nh=0.01）")
    ap.add_argument("--jobs", type=int, default=1, help="并行 worker 数")
    ap.add_argument("--track-progress", action="store_true", help="实时追踪任务进度，并写状态文件")

    ap.add_argument("tokens", nargs="*", help="如：1H0419-577 0148000201 0148000401  Mrk509 0130720101")

    args = ap.parse_args()

    out_root = args.out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    entries = scan_all_entries()
    full_idx = build_full_index(entries, recursive=args.recursive)
    index_txt = write_full_index_txt(out_root, entries, full_idx)
    print(f"[OK] wrote full index: {index_txt}")

    if not args.read_xcm:
        return

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
    obs_csv = dump_base / f"xspec_summary_obs_{src_tag}_{obs_tag}_{stamp}.csv"

    debug_dir = dump_base / "debug_xspec"
    if args.debug_xspec_out or args.track_progress:
        debug_dir.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    summary_rows_obs: List[Dict[str, Any]] = []

    lines.append("XSPEC param dump: keep RAW show par + PASS summary + ALL-param error table")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Selected targets: {len(selected)}")
    do_fit_effective = args.do_fit or args.do_error
    lines.append(f"do_fit={do_fit_effective}  do_error={args.do_error}  fit_iters={args.fit_iters}")
    lines.append(f"online_znh={'OFF' if args.no_online else 'ON'}")
    lines.append(f"jobs={args.jobs}")
    lines.append(f"track_progress={args.track_progress}")
    if args.do_error:
        lines.append(f"error_mode: {'FREE-only' if args.error_free_only else 'ALL(1-$npar)'}")
    lines.append("")

    znh_cache: Dict[str, Tuple[float, float, str, str]] = {}
    tasks: List[Dict[str, Any]] = []
    total_order = 0

    for e in selected:
        spec = e.path / "PN_spectrum_grp.fits"
        if not spec.exists():
            continue

        picks, note = full_idx.get((e.source, e.obsid), ([], ""))

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

        tom_slab = write_tom_xcm(e.path, e.source, e.obsid, slab=True, params=params, include_cpd=args.include_cpd)
        tom_sph = write_tom_xcm(e.path, e.source, e.obsid, slab=False, params=params, include_cpd=args.include_cpd)

        geom_to_pick: Dict[str, Optional[XcmPick]] = {"sphere": None, "slab": None}
        for p in sorted(picks, key=lambda x: x.mtime, reverse=True):
            if p.kind in geom_to_pick and geom_to_pick[p.kind] is None:
                geom_to_pick[p.kind] = p

        for geom in ("sphere", "slab"):
            pick = geom_to_pick[geom]
            if pick is None:
                summary_rows_obs.append(
                    build_empty_obs_row(
                        source=e.source,
                        obsid=e.obsid,
                        geom=geom,
                        z=z,
                        nh_1e22=nh_1e22,
                        simbad_used_name=sim_used_name,
                        pick_note=note,
                    )
                )
                continue

            rank = 1 if geom == "sphere" else 2
            tom = tom_sph if geom == "sphere" else tom_slab

            tasks.append({
                "order": total_order,
                "source": e.source,
                "obsid": e.obsid,
                "geom": geom,
                "rank": rank,
                "chi_kind": pick.kind,
                "obs_dir": str(e.path),
                "chi_path": str(pick.path),
                "tom_path": str(tom),
                "chi_mtime_str": pick.mtime_str,
                "pick_note": note,
                "simbad_used_name": sim_used_name,
                "z": z,
                "nh_1e22": nh_1e22,
                "w3nh_url": w3nh_url,
                "do_error": args.do_error,
                "do_fit": args.do_fit,
                "error_free_only": args.error_free_only,
                "fit_iters": args.fit_iters,
                "debug_xspec_out": args.debug_xspec_out,
                "track_progress": args.track_progress,
                "debug_dir": str(debug_dir) if (args.debug_xspec_out or args.track_progress) else "",
            })
            total_order += 1

    if args.track_progress:
        print(f"[INFO] total parallel tasks = {len(tasks)}")
        write_progress_summary(debug_dir, len(tasks), 0, 0, 0, [])

    results: List[Dict[str, Any]] = []
    jobs = max(1, int(args.jobs))

    finished_tags: List[str] = []
    failed_tasks = 0
    done_tasks = 0

    if jobs == 1:
        for i, t in enumerate(tasks, 1):
            if args.track_progress:
                print(f"[START {i}/{len(tasks)}] {t['source']} {t['obsid']} {t['geom']}")
            res = worker_run_fit_task(t)
            results.append(res)
            done_tasks += 1
            if not res["ok"]:
                failed_tasks += 1
            finished_tags.append(f"{res['source']} {res['obsid']} {res['geom']} ok={res['ok']}")
            if args.track_progress:
                print(f"[DONE  {done_tasks}/{len(tasks)}] {res['source']} {res['obsid']} {res['geom']} ok={res['ok']}")
                write_progress_summary(debug_dir, len(tasks), done_tasks, 0, failed_tasks, finished_tags)
    else:
        fut_map = {}
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for t in tasks:
                fut = ex.submit(worker_run_fit_task, t)
                fut_map[fut] = t
                if args.track_progress:
                    print(f"[SUBMIT] {t['source']} {t['obsid']} {t['geom']}")

            for fut in as_completed(fut_map):
                res = fut.result()
                results.append(res)
                done_tasks += 1
                if not res["ok"]:
                    failed_tasks += 1
                finished_tags.append(f"{res['source']} {res['obsid']} {res['geom']} ok={res['ok']}")
                if args.track_progress:
                    print(f"[DONE  {done_tasks}/{len(tasks)}] {res['source']} {res['obsid']} {res['geom']} ok={res['ok']}")
                    running_tasks = max(0, len(tasks) - done_tasks)
                    write_progress_summary(debug_dir, len(tasks), done_tasks, running_tasks, failed_tasks, finished_tags)

    results.sort(key=lambda x: x["order"])

    total_tables = 0
    for res in results:
        lines.extend(res["txt_lines"])
        summary_rows_obs.append(res["obs_row"])
        if res.get("has_table"):
            total_tables += 1

    dump_txt.write_text("\n".join(lines), encoding="utf-8")
    write_csv_rows(obs_csv, summary_rows_obs)

    print(f"[OK] wrote param dump: {dump_txt}")
    print(f"[OK] wrote summary obs csv: {obs_csv}")
    print(f"[INFO] parsed show-par tables: {total_tables}")
    print(f"[INFO] jobs used: {jobs} / system cpu_count={os.cpu_count()}")


if __name__ == "__main__":
    main()