#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import re
from datetime import datetime


# =============================================================================
# ✅ 留一个位置：默认搜寻文件夹（你只要改这一个日期目录就行）
# =============================================================================
DEFAULT_ROOT = Path("~/data/monk/plot/warmcorona/test/2026.02.02").expanduser().resolve()
DEFAULT_SEARCH_DIR = (DEFAULT_ROOT / "xmm_paramdump_with_error").resolve()
DEFAULT_OUT_DIR = (DEFAULT_ROOT / "latex_tables").resolve()

# =============================================================================
# ✅ 表格空隙：默认行距/列距（你嫌大/小就改这俩）
# =============================================================================
DEFAULT_ARRAYSTRETCH = 1.25   # 行距：1.15~1.35 常用
DEFAULT_TABCOLSEP = 8         # 列距(pt)：6~10 常用


SEP_RE = re.compile(r"^=+\s*$")
HDR_RE = re.compile(
    r"^(?P<src>[^/]+)/(?P<obsid>\d{10})\s+rank=(?P<rank>\d+)\s+chi_kind=(?P<kind>\w+)\s+chi_mtime=(?P<mtime>.+)$"
)

# 只读第二表：params(all)
PARAMS_ALL_RE = re.compile(r"^\s*params\(all\):", re.IGNORECASE)

# 第二表的行：par comp component parameter value pm errpair flags
# 注意 errpair 可能是 "--" 或 "-0.01,0.02" 或 "1,0.104101"（但后面会过滤掉 frozen/tied/redshift）
PARAMS_ROW_RE = re.compile(
    r"^\s*(?P<par>\d+)\s+(?P<comp>\d+)\s+(?P<component>\S+)\s+"
    r"(?P<parameter>.*?)\s+"
    r"(?P<value>[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"(?P<pm>--|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"(?P<errpair>--|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\s*,\s*[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
    r"(?P<flags>.*)$"
)

def safe_slug(s: str) -> str:
    out = []
    for ch in s:
        if ch.isalnum() or ch in "._-":
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_")


@dataclass
class ParRow:
    par: int
    comp: int
    component: str
    parameter: str
    value: float
    dminus: Optional[float]
    dplus: Optional[float]
    flags: str


@dataclass
class Block:
    src: str
    obsid: str
    rank: int
    kind: str
    mtime: str
    rows: List[ParRow]


# -------------------------
# 目录/文件搜寻
# -------------------------
def resolve_search_dir(root: Path, search_dir: Optional[Path]) -> Path:
    """
    优先使用 --search-dir；否则使用 root/xmm_paramdump_with_error
    """
    if search_dir is not None:
        return search_dir.expanduser().resolve()
    return (root.expanduser().resolve() / "xmm_paramdump_with_error").resolve()


def find_latest_dump(search_dir: Path) -> Path:
    search_dir = search_dir.expanduser().resolve()

    patterns = [
        "xmm_paramdump_with_error*.txt",
        "xspec_paramdump_*.txt",        # 兼容旧名字
        "*.txt",                        # 兜底
    ]

    cands: List[Path] = []
    for pat in patterns:
        cands.extend(search_dir.glob(pat))

    cands = [p for p in cands if p.is_file()]
    cands.sort(key=lambda p: p.stat().st_mtime, reverse=True)

    if not cands:
        raise FileNotFoundError(
            f"No dump txt found under: {search_dir}\n"
            f"Tried patterns: {patterns}"
        )
    return cands[0]


# -------------------------
# 解析：只读第二表 params(all)
# -------------------------
def parse_dump(txt: Path) -> List[Block]:
    lines = txt.read_text(encoding="utf-8", errors="ignore").splitlines()
    blocks: List[Block] = []
    i, n = 0, len(lines)

    while i < n:
        if not SEP_RE.match(lines[i]):
            i += 1
            continue
        if i + 1 >= n:
            break

        hm = HDR_RE.match(lines[i + 1].strip())
        if not hm:
            i += 1
            continue

        src = hm.group("src")
        obsid = hm.group("obsid")
        rank = int(hm.group("rank"))
        kind = hm.group("kind")
        mtime = hm.group("mtime").strip()

        j = i + 2

        # 找 params(all):
        while j < n and not PARAMS_ALL_RE.match(lines[j]):
            j += 1
        if j >= n:
            i = j
            continue

        # 找表头行 par comp（第二表）
        while j < n and not lines[j].lstrip().startswith("par comp"):
            j += 1
        if j >= n:
            i = j
            continue

        j += 1
        rows: List[ParRow] = []

        while j < n:
            s = lines[j]
            if not s.strip():
                break
            if SEP_RE.match(s):
                break

            rm = PARAMS_ROW_RE.match(s)
            if rm:
                par_i = int(rm.group("par"))
                comp_i = int(rm.group("comp"))
                component = rm.group("component").strip()
                parameter = rm.group("parameter").strip()
                value = float(rm.group("value"))
                errpair = rm.group("errpair").strip()
                flags = rm.group("flags").strip()

                dminus = None
                dplus = None

                # ✅ 只认 err:4float(d-,d+,lo,hi) 的这一对（就是你要的第二表 error）
                # 其它 err:2float(...) 或 frozen/tied 都会被 keep_row 过滤掉
                if "err:4float" in flags.lower() and errpair != "--" and "," in errpair:
                    a, b = errpair.split(",", 1)
                    try:
                        dminus = float(a.strip())
                        dplus = float(b.strip())
                    except Exception:
                        dminus, dplus = None, None

                rows.append(ParRow(
                    par=par_i, comp=comp_i,
                    component=component, parameter=parameter,
                    value=value, dminus=dminus, dplus=dplus,
                    flags=flags
                ))
            j += 1

        blocks.append(Block(src=src, obsid=obsid, rank=rank, kind=kind, mtime=mtime, rows=rows))
        i = j

    return blocks


# -------------------------
# 模型/参数名 & 过滤
# -------------------------
def model_name_from_component(comp: str) -> str:
    c = comp.lower()
    if c.startswith("warmcom_"):
        return "warmcom"
    if c.startswith("zgauss"):
        return "zgauss"
    if c.startswith("ztbabs"):
        return "zTBabs"
    if c.startswith("tbabs"):
        return "TBabs"
    if c.startswith("nthcomp"):
        return "nthComp"
    if c.startswith("xillver"):
        return "xillver"
    return comp

def is_tied_or_frozen(flags: str) -> bool:
    low = flags.lower()
    if "frozen" in low:
        return True
    # 你的 flags 里 tied 是 "= p2" / "= p4" 这种
    if "= p" in low or "=p" in low:
        return True
    return False

def is_redshift_like(comp: str, param: str) -> bool:
    c = comp.lower()
    p = param.lower().strip()
    if "redshift" in p:
        return True
    if c.startswith("warmcom_") and c.endswith("z"):
        return True
    if c == "xillver" and p == "z":
        return True
    if c.startswith("zgauss") and p == "redshift":
        return True
    return False

def keep_row(row: ParRow) -> bool:
    # ✅ 你不想要的：redshift 相关 & frozen/tied
    if is_redshift_like(row.component, row.parameter):
        return False
    if is_tied_or_frozen(row.flags):
        return False
    return True

def pretty_param_name(row: ParRow) -> str:
    c = row.component.lower()
    p = row.parameter.strip()
    pl = p.lower()

    if c.startswith("warmcom_"):
        if c.endswith("te"):
            return r"{\small $kT_{\rm e}$ (keV)}"
        if c.endswith("tau"):
            return r"{\small $\tau$}"
        if c.endswith("norm"):
            return r"{\small norm}"
        return r"{\small " + row.component.replace("_", r"\_") + r"}"

    if pl == "gamma":
        return r"{\small $\Gamma$}"
    if pl.startswith("linee"):
        return r"{\small LineE (keV)}"
    if pl == "sigma":
        return r"{\small $\sigma$ (keV)}"
    if pl.startswith("nh"):
        if row.component.lower().startswith("tbabs"):
            return r"{\small $N_{\rm H,Gal}$ [$10^{22}$]}"
        if row.component.lower().startswith("ztbabs"):
            return r"{\small $N_{\rm H,int}$ [$10^{22}$]}"
        return r"{\small $N_{\rm H}$ [$10^{22}$]}"
    if pl == "afe":
        return r"{\small $A_{\rm Fe}$}"
    if pl == "logxi":
        return r"{\small $\log\xi$}"
    if pl.startswith("ecut"):
        return r"{\small $E_{\rm cut}$ (keV)}"
    if pl.startswith("incl"):
        return r"{\small Incl (deg)}"
    if pl == "refl_frac":
        return r"{\small refl\_frac}"
    if pl == "norm":
        return r"{\small norm}"

    return r"{\small " + p.replace("_", r"\_") + r"}"


# -------------------------
# 数值格式：只用 d-/d+（第二表），并确保数学模式
# -------------------------
def fmt_num(x: float, sig: int = 4) -> str:
    ax = abs(x)
    if ax == 0:
        return "0"
    if ax < 1e-3 or ax >= 1e4:
        return f"{x:.{sig}e}".replace("e+0", "e+").replace("e-0", "e-")
    return f"{x:.{sig}g}"

def ensure_math(s: str) -> str:
    # ✅ 避免 ^ _ 报错：强制进数学模式
    return r"\ensuremath{" + s + r"}"

def build_cell(row: Optional[ParRow]) -> str:
    if row is None:
        return r"--"

    # ✅ 有 error（d-/d+）就用它；否则只给中心值
    if row.dminus is not None and row.dplus is not None:
        v = fmt_num(row.value)
        dm = fmt_num(row.dminus)
        dp = fmt_num(row.dplus)

        if not dp.startswith(("+", "-")):
            dp = "+" + dp
        if not dm.startswith(("+", "-")):
            dm = "-" + dm

        return ensure_math(rf"{v}^{{{dp}}}_{{{dm}}}")

    return ensure_math(fmt_num(row.value))


# -------------------------
# 渲染 LaTeX
# -------------------------
def render_one_source(blocks: List[Block], source: str, kinds: List[str],
                      arraystretch: float, tabcolsep: int) -> str:
    bsel = [b for b in blocks if b.src.lower() == source.lower()]
    if not bsel:
        raise RuntimeError(f"No blocks for source={source}")

    kind_order = {k.lower(): i for i, k in enumerate(kinds)}
    bsel.sort(key=lambda b: (b.obsid, kind_order.get(b.kind.lower(), 99), b.rank))

    col_titles = [f"{b.obsid}-{b.kind}" for b in bsel]
    colspec = "ll" + "c" * len(col_titles)

    # ✅ warmcom 归一：te/tau/norm 只出一行
    def row_key(r: ParRow) -> Tuple[str, str]:
        c = r.component.lower()
        if c.startswith("warmcom_"):
            if c.endswith("te"):
                return ("warmcom", "te")
            if c.endswith("tau"):
                return ("warmcom", "tau")
            if c.endswith("norm"):
                return ("warmcom", "norm")
            return ("warmcom", c)
        return (r.component, " ".join(r.parameter.lower().split()))

    all_rows: List[ParRow] = []
    for b in bsel:
        for r in b.rows:
            if keep_row(r):
                all_rows.append(r)

    all_rows.sort(key=lambda r: (model_name_from_component(r.component).lower(), r.par))

    row_defs: List[Tuple[str, str, Tuple[str, str]]] = []
    seen = set()
    for r in all_rows:
        model = model_name_from_component(r.component)
        pname = pretty_param_name(r)
        key = row_key(r)
        sig = (model, pname, key)
        if sig in seen:
            continue
        seen.add(sig)
        row_defs.append((model, pname, key))

    lines: List[str] = []
    lines.append(r"\begin{table*}[t]")
    lines.append(r"\centering")
    # ✅ 空隙更大
    lines.append(rf"\renewcommand{{\arraystretch}}{{{arraystretch}}}")
    lines.append(rf"\setlength{{\tabcolsep}}{{{tabcolsep}pt}}")
    lines.append(rf"\caption{{Best-fit parameters for {source}.}}")
    lines.append(rf"\label{{tab:{safe_slug(source)}}}")
    lines.append(rf"\begin{{tabular}}{{{colspec}}}")
    lines.append(r"\hline")
    lines.append(r"Model & Parameter & " + " & ".join(col_titles) + r" \\")
    lines.append(r"\hline")

    last_model = None
    for model, pname, key in row_defs:
        show_model = model if model != last_model else ""
        last_model = model

        row_vals: List[str] = []
        for b in bsel:
            hit = None
            for rr in b.rows:
                if not keep_row(rr):
                    continue
                if row_key(rr) == key:
                    hit = rr
                    break
            row_vals.append(build_cell(hit))

        lines.append(f"{show_model} & {pname} & " + " & ".join(row_vals) + r" \\")
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table*}")
    lines.append("")
    return "\n".join(lines)


def build_out_path(out_dir: Path, source: str, obsids: List[str], base_tag: str) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    obsids_sorted = sorted(set(obsids))
    if len(obsids_sorted) == 1:
        obs_tag = obsids_sorted[0]
    else:
        obs_tag = f"{obsids_sorted[0]}-{obsids_sorted[-1]}_{len(obsids_sorted)}obs"
    return out_dir / f"table_{safe_slug(source)}_{obs_tag}_{base_tag}_{stamp}.tex"


def main():
    ap = argparse.ArgumentParser(
        description="Parse XSPEC dump txt and generate LaTeX table using ONLY params(all) error (2nd table)."
    )
    ap.add_argument("--txt", type=Path, default=None, help="dump txt path")
    ap.add_argument("--auto-latest", action="store_true", help="auto find latest dump txt")

    # ✅ root 默认就是你留的位置
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help="root test dir (default is the hard-coded DEFAULT_ROOT in script)")

    ap.add_argument("--search-dir", type=Path, default=None,
                    help="override search dir (default: <root>/xmm_paramdump_with_error)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="override output dir (default: <root>/latex_tables)")

    ap.add_argument("--kinds", nargs="*", default=["sphere", "slab"], help="column order by kind")
    ap.add_argument("--arraystretch", type=float, default=DEFAULT_ARRAYSTRETCH, help="LaTeX row spacing")
    ap.add_argument("--tabcolsep", type=int, default=DEFAULT_TABCOLSEP, help="LaTeX column spacing (pt)")

    args = ap.parse_args()

    root = args.root.expanduser().resolve()
    search_dir = resolve_search_dir(root, args.search_dir)

    out_dir = (args.out_dir.expanduser().resolve() if args.out_dir else (root / "latex_tables").resolve())
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.auto_latest:
        txt_path = find_latest_dump(search_dir)
    else:
        if args.txt is None:
            raise SystemExit("Please provide --txt, or use --auto-latest (optionally with --root/--search-dir).")
        txt_path = args.txt.expanduser().resolve()

    blocks = parse_dump(txt_path)
    if not blocks:
        raise SystemExit(f"[ERROR] No blocks parsed from: {txt_path}")

    base_tag = safe_slug(txt_path.stem)
    if not base_tag:
        base_tag = "dump"

    bysrc: Dict[str, List[Block]] = {}
    for b in blocks:
        bysrc.setdefault(b.src, []).append(b)

    for src, blist in sorted(bysrc.items(), key=lambda kv: kv[0].lower()):
        tex = render_one_source(
            blist, src, kinds=args.kinds,
            arraystretch=args.arraystretch,
            tabcolsep=args.tabcolsep
        )
        out_path = build_out_path(out_dir, src, [bb.obsid for bb in blist], base_tag=base_tag)
        out_path.write_text(tex, encoding="utf-8")
        print(f"[OK] wrote: {out_path}")

    print(f"[INFO] from dump : {txt_path}")
    print(f"[INFO] search_dir: {search_dir}")
    print(f"[INFO] out_dir   : {out_dir}")


if __name__ == "__main__":
    main()
