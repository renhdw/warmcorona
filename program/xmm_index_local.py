#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import re
from datetime import datetime
import os


# =========================
# 0) 环境路径（按需改）
# =========================
XMM_BASE_DIR = Path("/home/hdw/data/XMM").resolve()
DEFAULT_OUT_ROOT = Path("~/data/monk/plot/warmcorona/test/2026.03.24/xcm_upgrade").expanduser().resolve()

# ---------- warmcom sphere ----------
SPHERE_MODEL_DIR  = os.path.expanduser("~/data/monk/plot/warmcorona/model/warmcom/warmcom_sphere_log/warmcom_sphere")
SPHERE_PKG        = "warmcom_sphere"
SPHERE_MODEL_NAME = "warmcomsphere"

# ---------- warmcom slab ----------
SLAB_MODEL_DIR  = os.path.expanduser("~/data/monk/plot/warmcorona/model/warmcom/warmcom_slab_log/warmcom_slab")
SLAB_PKG        = "warmcom_slab"
SLAB_MODEL_NAME = "warmcomslab"

# 输出命名
CALL_SLAB_TEMPLATE   = "fit_model_local_slab_{source}_{obsid}.xcm"
CALL_SPHERE_TEMPLATE = "fit_model_local_sphere_{source}_{obsid}.xcm"
CHI_LOCAL_SUFFIX = "_local"

_EXCLUDE_CHI_SUBSTRS = ["gdiskwien"]

_OBSID_RE = re.compile(r"^\d{10}$")
_CHI_XCM_RE = re.compile(r"^chi_.*\.xcm$", re.IGNORECASE)

# 参数行：既支持简单 call xcm，也支持 chi xcm 里 "= p2" 这种行
_PARAM_LINE_RE = re.compile(r"^\s*(?:[-+]?[\d.]+|=\s*p?\d+)", re.IGNORECASE)
_INT_RE = re.compile(r"^\d+$")
_RANGE_RE = re.compile(r"^(\d+)-(\d+)$")


# =========================
# 1) 数据结构
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
    kind: str  # slab / sphere

    @property
    def mtime_str(self) -> str:
        return datetime.fromtimestamp(self.mtime).strftime("%Y-%m-%d %H:%M:%S")


# =========================
# 2) 小工具
# =========================
def _is_source_dir(p: Path) -> bool:
    return p.is_dir() and not p.name.startswith(".")


def _is_obsid_dir(p: Path) -> bool:
    return p.is_dir() and (_OBSID_RE.match(p.name) is not None)


def should_exclude_chi(path: Path) -> bool:
    low = path.name.lower()
    return any(sub in low for sub in _EXCLUDE_CHI_SUBSTRS)


def read_text_safe(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def guess_kind_from_name_or_text(name: str, text: str = "") -> str:
    low_name = name.lower()
    low_text = text.lower()

    if "slab" in low_name:
        return "slab"
    if "sphere" in low_name:
        return "sphere"

    if "warmcomslab" in low_text or "warmcom_slab" in low_text:
        return "slab"
    if "warmcomsphere" in low_text or "warmcom_sphere" in low_text:
        return "sphere"

    if "pure_tom_slab" in low_text or "warmcom_0.1-0.6_5-25" in low_text:
        return "slab"
    if "tom_new_107" in low_text or "smoothed_tv_tom_new_107" in low_text:
        return "sphere"

    return "unknown"


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


def expand_specs(spec_tokens: List[str]) -> List[int]:
    nums: List[int] = []
    for tok in spec_tokens:
        tok = tok.strip()
        if not tok:
            continue
        m = _RANGE_RE.match(tok)
        if m:
            a = int(m.group(1))
            b = int(m.group(2))
            if a <= b:
                nums.extend(range(a, b + 1))
            else:
                nums.extend(range(a, b - 1, -1))
            continue
        if _INT_RE.match(tok):
            nums.append(int(tok))
    return nums


# =========================
# 3) 参数编号映射
# =========================
def remap_par_num_case_simple(n: int) -> int:
    """
    旧:
      1 zTBabs.nH
      2 zTBabs.z
      3 TBabs.nH
      4...
    新:
      1 TBabs.nH
      2 zTBabs.nH
      3 zTBabs.z
      4...
    """
    if n == 1:
        return 2
    if n == 2:
        return 3
    if n == 3:
        return 1
    return n


def remap_par_num_case_wa(n: int) -> int:
    """
    旧:
      1 zTBabs.nH
      2 zTBabs.z
      3 wa.PARAM1
      4 wa.PARAM2
      5 wa.z
      6 TBabs.nH
      7...
    新:
      1 TBabs.nH
      2 zTBabs.nH
      3 zTBabs.z
      4 wa.PARAM1
      5 wa.PARAM2
      6 wa.z
      7...
    """
    if n == 6:
        return 1
    if n == 1:
        return 2
    if n == 2:
        return 3
    if n == 3:
        return 4
    if n == 4:
        return 5
    if n == 5:
        return 6
    return n


def remap_par_num_by_case(n: int, case: str) -> int:
    if case == "simple":
        return remap_par_num_case_simple(n)
    if case == "wa":
        return remap_par_num_case_wa(n)
    return n


def remap_link_target_only(text: str, remap_func) -> str:
    def repl_eq_num(m):
        return "=" + str(remap_func(int(m.group(1))))

    def repl_eq_num_sp(m):
        return "= " + str(remap_func(int(m.group(1))))

    def repl_eq_p(m):
        return "=p" + str(remap_func(int(m.group(1))))

    def repl_eq_p_sp(m):
        return "= p" + str(remap_func(int(m.group(1))))

    text = re.sub(r"=(\d+)\b", repl_eq_num, text)
    text = re.sub(r"=\s+(\d+)\b", repl_eq_num_sp, text)
    text = re.sub(r"=p(\d+)\b", repl_eq_p, text, flags=re.I)
    text = re.sub(r"=\s+p(\d+)\b", repl_eq_p_sp, text, flags=re.I)
    return text


# =========================
# 4) 扫描 obs
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
# 5) 只取最新 chi_slab / chi_sphere
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
            xcm_dir = obs_dir / "xcm"
            if xcm_dir.exists() and xcm_dir.is_dir():
                for p in xcm_dir.iterdir():
                    if p.is_file() and _CHI_XCM_RE.match(p.name) and (not should_exclude_chi(p)):
                        out.append(p)

    uniq: Dict[Path, Path] = {}
    for p in out:
        uniq[p.resolve()] = p
    return list(uniq.values())


def pick_latest_chi_slab_sphere(obs_dir: Path, recursive: bool) -> Tuple[List[XcmPick], str]:
    allp: List[XcmPick] = []

    for f in list_chi_xcms(obs_dir, recursive):
        try:
            st = f.stat()
            txt = read_text_safe(f)
        except Exception:
            continue

        kind = guess_kind_from_name_or_text(f.name, txt)
        if kind not in {"slab", "sphere"}:
            continue

        allp.append(XcmPick(path=f, mtime=st.st_mtime, kind=kind))

    if not allp:
        return [], "NO valid chi_*.xcm found"

    slab = sorted([p for p in allp if p.kind == "slab"], key=lambda x: x.mtime, reverse=True)
    sphere = sorted([p for p in allp if p.kind == "sphere"], key=lambda x: x.mtime, reverse=True)

    picks: List[XcmPick] = []
    note_parts: List[str] = []

    if sphere:
        picks.append(sphere[0])
        note_parts.append("latest sphere picked")
    else:
        note_parts.append("no sphere")

    if slab:
        picks.append(slab[0])
        note_parts.append("latest slab picked")
    else:
        note_parts.append("no slab")

    picks.sort(key=lambda x: x.mtime, reverse=True)
    return picks, "; ".join(note_parts)


def build_full_index(entries: List[ObsEntry], recursive: bool) -> Dict[Tuple[str, str], Tuple[List[XcmPick], str]]:
    idx: Dict[Tuple[str, str], Tuple[List[XcmPick], str]] = {}
    for e in entries:
        idx[(e.source, e.obsid)] = pick_latest_chi_slab_sphere(e.path, recursive)
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
    lines.append("XMM FULL INDEX (latest chi slab/sphere only)")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"XMM_BASE_DIR: {XMM_BASE_DIR}")
    lines.append("Rule: only latest slab chi + latest sphere chi")
    lines.append("")

    for src, obslist in bysrc.items():
        lines.append(f"== {src} (obsids={len(obslist)}) ==")
        for e in obslist:
            picks, note = idx.get((e.source, e.obsid), ([], ""))
            lines.append(f"{src}/{e.obsid}  picked={len(picks)}  note={note}")
            if not picks:
                lines.append("  (NO valid chi_*.xcm found)")
            else:
                for k, p in enumerate(sorted(picks, key=lambda x: x.mtime, reverse=True), 1):
                    lines.append(f"  {k}) kind={p.kind:<7s}  {p.path}  mtime={p.mtime_str}")
        lines.append("")

    out_txt.write_text("\n".join(lines), encoding="utf-8")
    return out_txt


# =========================
# 6) 找对应的调用 xcm
# =========================
def list_call_xcms_in_xcm_dir(obs_dir: Path) -> List[Path]:
    xcm_dir = obs_dir / "xcm"
    if not xcm_dir.exists() or not xcm_dir.is_dir():
        return []

    out: List[Path] = []
    for p in xcm_dir.iterdir():
        if not p.is_file():
            continue
        if p.suffix.lower() != ".xcm":
            continue
        if p.name.lower().startswith("chi_"):
            continue
        out.append(p)
    return out


def pick_latest_call_xcm(obs_dir: Path, kind: str) -> Optional[Path]:
    cands = []
    for p in list_call_xcms_in_xcm_dir(obs_dir):
        txt = read_text_safe(p)
        k = guess_kind_from_name_or_text(p.name, txt)
        if k == kind:
            cands.append(p)

    if not cands:
        return None

    cands.sort(key=lambda x: x.stat().st_mtime, reverse=True)
    return cands[0]


# =========================
# 7) model / 参数块识别
# =========================
def split_model_block(lines: List[str]) -> Tuple[List[str], Optional[int], Optional[int], Optional[int]]:
    model_idx = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("model "):
            model_idx = i
            break

    if model_idx is None:
        return lines[:], None, None, None

    j = model_idx + 1
    while j < len(lines) and lines[j].strip() == "":
        j += 1
    param_start = j

    while j < len(lines):
        s = lines[j].strip()
        if s == "":
            j += 1
            continue
        if _PARAM_LINE_RE.match(s):
            j += 1
            continue
        break

    return lines[:model_idx], model_idx, param_start, j


def detect_model_case(old_model_line: str) -> str:
    s = old_model_line.lower().replace(" ", "")

    if "ztbabs*mtable{" in s and "}*tbabs(" in s:
        return "wa"

    if "ztbabs*tbabs(" in s:
        return "simple"

    return "other"


def build_local_lmod_line(kind: str) -> str:
    if kind == "sphere":
        return f"lmod {SPHERE_PKG} {SPHERE_MODEL_DIR}"
    if kind == "slab":
        return f"lmod {SLAB_PKG} {SLAB_MODEL_DIR}"
    raise ValueError(f"Unknown kind: {kind}")


def strip_old_local_lmod(pre_lines: List[str]) -> List[str]:
    out = []
    for line in pre_lines:
        low = line.strip().lower()
        if low.startswith("lmod ") and ("warmcom_sphere" in low or "warmcom_slab" in low):
            continue
        out.append(line)
    return out


def ensure_local_lmod(pre_lines: List[str], kind: str) -> List[str]:
    out = strip_old_local_lmod(pre_lines)
    want = build_local_lmod_line(kind)

    for x in out:
        if x.strip() == want:
            return out

    insert_at = None
    for i, line in enumerate(out):
        if line.strip().lower().startswith("lmod relxill"):
            insert_at = i + 1
            break

    if insert_at is None:
        out.append(want)
    else:
        out.insert(insert_at, want)
    return out


def build_local_model_line_from_old(old_model_line: str, kind: str) -> str:
    warm_name = SLAB_MODEL_NAME if kind == "slab" else SPHERE_MODEL_NAME
    s = old_model_line.strip()
    s = re.sub(r"^\s*model\s+", "", s, flags=re.I)

    # atable -> local model
    s = re.sub(r"atable\s*\{[^}]+\}", warm_name, s, flags=re.I)

    case = detect_model_case(old_model_line)

    if case == "simple":
        # zTBabs*TBabs(...) -> TBabs*zTBabs(...)
        s = re.sub(
            r"^\s*zTBabs\s*\*\s*TBabs\s*\(",
            r"TBabs*zTBabs(",
            s,
            flags=re.I
        )
    elif case == "wa":
        # zTBabs*mtable{...}*TBabs(...) -> TBabs*zTBabs*mtable{...}(...)
        s = re.sub(
            r"^\s*zTBabs\s*\*\s*(mtable\s*\{[^}]+\})\s*\*\s*TBabs\s*\(",
            r"TBabs*zTBabs*\1(",
            s,
            flags=re.I
        )

    return "model " + s


def reorder_param_lines(old_param_lines: List[str], case: str) -> List[str]:
    core = [x for x in old_param_lines if x.strip() != ""]

    if case == "simple":
        if len(core) < 3:
            return old_param_lines[:]
        # old [1,2,3,...] -> new [3,1,2,...]
        return [core[2], core[0], core[1]] + core[3:]

    if case == "wa":
        if len(core) < 6:
            return old_param_lines[:]
        # old [1,2,3,4,5,6,...] -> new [6,1,2,3,4,5,...]
        return [core[5], core[0], core[1], core[2], core[3], core[4]] + core[6:]

    return old_param_lines[:]


def remap_param_block_lines(param_lines: List[str], case: str) -> List[str]:
    """
    参数块内部也可能有纯 link 行，比如:
        = p2
        =p2
    这些必须一起 remap 成新的参数号。
    """
    out = []
    for line in param_lines:
        s = line.rstrip("\n")
        s = remap_link_target_only(s, lambda x: remap_par_num_by_case(x, case))
        out.append(s)
    return out


# =========================
# 8) 后处理命令 remap
# =========================
def remap_freeze_thaw_line(line: str, case: str) -> str:
    stripped = line.strip()
    if not stripped:
        return line.rstrip("\n")

    parts = stripped.split()
    cmd = parts[0]
    specs = parts[1:]
    nums = expand_specs(specs)
    if not nums:
        return line.rstrip("\n")

    nums_new = [remap_par_num_by_case(x, case) for x in nums]
    return f"{cmd} {compress_ranges(nums_new)}"


def remap_new_line(line: str, case: str) -> str:
    s = line.rstrip("\n")
    m = re.match(r"^(\s*new\s+)(\d+)(.*)$", s, flags=re.I)
    if not m:
        return s

    prefix = m.group(1)
    pnum = int(m.group(2))
    tail = m.group(3)

    pnum_new = remap_par_num_by_case(pnum, case) if pnum != 0 else 0
    tail = remap_link_target_only(tail, lambda x: remap_par_num_by_case(x, case))
    return f"{prefix}{pnum_new}{tail}"


def remap_tie_line(line: str, case: str) -> str:
    s = remap_link_target_only(line.rstrip("\n"), lambda x: remap_par_num_by_case(x, case))
    m = re.match(r"^(\s*)(\d+)(\s*=.*)$", s)
    if m:
        p = int(m.group(2))
        return f"{m.group(1)}{remap_par_num_by_case(p, case)}{m.group(3)}"
    return s


def remap_post_model_line(line: str, case: str) -> str:
    low = line.strip().lower()
    if low.startswith("freeze ") or low.startswith("thaw ") or low.startswith("untie "):
        return remap_freeze_thaw_line(line, case)
    if low.startswith("new "):
        return remap_new_line(line, case)
    if low.startswith("tie ") or re.match(r"^\s*\d+\s*=", line):
        return remap_tie_line(line, case)
    return line.rstrip("\n")


# =========================
# 9) 生成新 xcm
# =========================
def migrate_call_xcm_text(old_text: str, kind: str) -> str:
    """
    call xcm 只改：
      - 加 local lmod
      - 改 model 行
    其余保持不动
    """
    lines = old_text.splitlines()
    pre, model_idx, _, _ = split_model_block(lines)

    if model_idx is None:
        pre2 = ensure_local_lmod(lines[:], kind)
        return "\n".join(pre2).rstrip() + "\n"

    out = lines[:]
    out_pre = ensure_local_lmod(out[:model_idx], kind)
    out_model = build_local_model_line_from_old(out[model_idx], kind)

    final = []
    final.extend(out_pre)
    final.append(out_model)
    final.extend(out[model_idx + 1:])
    return "\n".join(final).rstrip() + "\n"


def migrate_full_chi_text(old_text: str, kind: str) -> str:
    """
    chi_xcm:
      - 加 local lmod
      - 改 model
      - 调整参数块顺序
      - remap 参数块内部 link 行
      - remap freeze/new/tie
    """
    lines = old_text.splitlines()
    pre, model_idx, param_start, param_end = split_model_block(lines)

    if model_idx is None or param_start is None or param_end is None:
        pre2 = ensure_local_lmod(lines[:], kind)
        return "\n".join(pre2).rstrip() + "\n"

    old_model_line = lines[model_idx]
    case = detect_model_case(old_model_line)

    old_param_lines = lines[param_start:param_end]
    post_lines = lines[param_end:]

    pre2 = ensure_local_lmod(pre, kind)
    new_model_line = build_local_model_line_from_old(old_model_line, kind)
    new_param_lines = reorder_param_lines(old_param_lines, case)
    new_param_lines = remap_param_block_lines(new_param_lines, case)
    post2 = [remap_post_model_line(x, case) for x in post_lines]

    out = []
    out.extend(pre2)
    out.append(new_model_line)
    out.extend(new_param_lines)
    out.extend(post2)
    return "\n".join(out).rstrip() + "\n"


def output_call_name(source: str, obsid: str, kind: str) -> str:
    if kind == "sphere":
        return CALL_SPHERE_TEMPLATE.format(source=source, obsid=obsid)
    return CALL_SLAB_TEMPLATE.format(source=source, obsid=obsid)


def output_chi_local_name(old_path: Path) -> str:
    return old_path.stem + CHI_LOCAL_SUFFIX + old_path.suffix


# =========================
# 10) targets
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
# 11) main
# =========================
def main():
    ap = argparse.ArgumentParser(
        description=(
            "Scan XMM dir, always build index, pick latest chi_slab/chi_sphere, "
            "upgrade model-call xcm and chi_xcm to localmodel versions. No XSPEC run."
        )
    )

    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    ap.add_argument("--recursive", action="store_true", help="递归搜 chi_*.xcm")

    ap.add_argument("--all", action="store_true", help="对所有 source/obsid 处理")
    ap.add_argument("tokens", nargs="*", help="如：RE1034+396 0675440301")

    ap.add_argument("--make-call-xcm", action="store_true", help="生成 local 模型调用 xcm")
    ap.add_argument("--make-chi-xcm", action="store_true", help="生成 local 版 chi_xcm")
    ap.add_argument("--write-into-xcm-dir", action="store_true", help="call_xcm 写入 obsid/xcm/；chi_xcm 写回原 chi 同级目录")
    ap.add_argument("--overwrite", action="store_true", help="允许覆盖已有文件")

    args = ap.parse_args()

    out_root = args.out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    entries = scan_all_entries()
    full_idx = build_full_index(entries, recursive=args.recursive)
    index_txt = write_full_index_txt(out_root, entries, full_idx)
    print(f"[OK] wrote full index: {index_txt}")

    if not args.make_call_xcm and not args.make_chi_xcm:
        print("[INFO] no generation flag given. index only.")
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

    if not selected:
        print("[WARN] no targets selected.")
        return

    log_lines: List[str] = []
    log_lines.append("XCM UPGRADE LOG")
    log_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log_lines.append(f"make_call_xcm={args.make_call_xcm}  make_chi_xcm={args.make_chi_xcm}")
    log_lines.append(f"write_into_xcm_dir={args.write_into_xcm_dir}  overwrite={args.overwrite}")
    log_lines.append("")

    n_call = 0
    n_chi = 0

    for e in selected:
        picks, note = full_idx.get((e.source, e.obsid), ([], ""))
        if not picks:
            log_lines.append(f"{e.source}/{e.obsid}  SKIP  no valid chi_xcm   note={note}")
            continue

        log_lines.append("=" * 100)
        log_lines.append(f"{e.source}/{e.obsid}  picked={len(picks)}  note={note}")

        for rank, pick in enumerate(sorted(picks, key=lambda x: x.mtime, reverse=True), 1):
            try:
                chi_old_text = read_text_safe(pick.path)
            except Exception as ex:
                log_lines.append(f"  rank={rank}  READ CHI FAIL  {pick.path}  err={ex}")
                continue

            kind = pick.kind
            log_lines.append(f"  rank={rank}  kind={kind}  chi={pick.path}")

            # ---------- call xcm ----------
            if args.make_call_xcm:
                call_old = pick_latest_call_xcm(e.path, kind)
                if call_old is None:
                    log_lines.append(f"    [SKIP] no call xcm found for kind={kind} in {e.path / 'xcm'}")
                else:
                    try:
                        call_old_text = read_text_safe(call_old)
                        call_new_text = migrate_call_xcm_text(call_old_text, kind)
                    except Exception as ex:
                        log_lines.append(f"    [FAIL] call migrate failed: {call_old} err={ex}")
                        call_new_text = None

                    if call_new_text is not None:
                        if args.write_into_xcm_dir:
                            call_out_dir = e.path / "xcm"
                        else:
                            call_out_dir = out_root / "upgraded_xcm" / e.source / e.obsid / "xcm"
                        call_out_dir.mkdir(parents=True, exist_ok=True)

                        call_name = output_call_name(e.source, e.obsid, kind)
                        call_out = call_out_dir / call_name

                        if call_out.exists() and (not args.overwrite):
                            log_lines.append(f"    [SKIP] call exists: {call_out}")
                        else:
                            call_out.write_text(call_new_text, encoding="utf-8")
                            n_call += 1
                            log_lines.append(f"    [OK] call -> {call_out}  (from {call_old.name})")

            # ---------- chi xcm ----------
            if args.make_chi_xcm:
                try:
                    chi_new_text = migrate_full_chi_text(chi_old_text, kind)
                except Exception as ex:
                    log_lines.append(f"    [FAIL] chi migrate failed: {pick.path} err={ex}")
                    chi_new_text = None

                if chi_new_text is not None:
                    if args.write_into_xcm_dir:
                        chi_out_dir = pick.path.parent
                    else:
                        chi_out_dir = out_root / "upgraded_xcm" / e.source / e.obsid / "chi"
                    chi_out_dir.mkdir(parents=True, exist_ok=True)

                    chi_name = output_chi_local_name(pick.path)
                    chi_out = chi_out_dir / chi_name

                    if chi_out.exists() and (not args.overwrite):
                        log_lines.append(f"    [SKIP] chi exists: {chi_out}")
                    else:
                        chi_out.write_text(chi_new_text, encoding="utf-8")
                        n_chi += 1
                        log_lines.append(f"    [OK] chi  -> {chi_out}")

    log_dir = out_root / "upgrade_log"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_txt = log_dir / f"xcm_upgrade_log_{stamp}.txt"

    log_lines.append("")
    log_lines.append(f"TOTAL call_xcm written: {n_call}")
    log_lines.append(f"TOTAL chi_xcm  written: {n_chi}")
    log_txt.write_text("\n".join(log_lines), encoding="utf-8")

    print(f"[OK] wrote upgrade log: {log_txt}")
    print(f"[INFO] call_xcm written: {n_call}")
    print(f"[INFO] chi_xcm  written: {n_chi}")


if __name__ == "__main__":
    main()