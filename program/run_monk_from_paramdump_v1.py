#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_monk_from_paramdump_v1.py  (fixed & stable)

功能（按你现在的工作流）：
1) 在 DEFAULT_SEARCH_DIR 里找最新的 xspec_paramdump_with_error_*.txt
2) 解析一个 txt 里可能包含的多个 SRC/OBSID 段落（每段带 geom=slab|sphere）
3) 从每段里抓 warmcom 拟合参数：te, tau, z, norm
4) 为每个 task 创建固定目录结构（你满意的那种）：
   DATA_DIR/<SRC>/<OBSID>/te_xxx_tau_xxx/
       ├─ sphere/   或 slab/          (只会有一个几何)
       ├─ calspec/
       ├─ logs/
       ├─ record.json
       ├─ xspec_warmcom_model.xcm
       ├─ xspec_warmcom_model.qdp
       └─ compare_warmcom_vs_monk.pdf
5) 跑 MONK：
   - 如果 calspec/en.dat + flux.dat 已存在且非空 => 默认跳过（不覆盖你已有数据）
   - 否则：复制模板 params -> (sphere|slab)/params.txt，写 te/tau，然后跑 sphere/slab，再跑 calspec
6) 跑 XSPEC 导出 warmcom-only 的 qdp：
   - 使用 `model atable{...}` + 4 行数值输入（不使用 newpar）
   - 加 show par 做 sanity check，确保参数真的设上了
7) 画图：上面双谱（XSPEC vs MONK），下面 ratio（XSPEC/MONK）
   ✅ 本版本：MONK 原始谱会乘上 paramdump 里的 norm 后再对比（与你 XSPEC QDP 对齐）

用法（保持简单）：
  python3 run_monk_from_paramdump_v1.py
  python3 run_monk_from_paramdump_v1.py --debug-parse
  python3 run_monk_from_paramdump_v1.py --only-geom sphere
  python3 run_monk_from_paramdump_v1.py --paramdump /path/to/xxx.txt

v1-binary+norm-fix:
- calspec/en.dat & flux.dat 支持 binary double(<d) 读取（优先按二进制读，失败再按文本读）
- XSPEC export 使用 bytes capture + replace decode，避免任何 UnicodeDecodeError
- Plot 时：MONK flux *= norm（对齐 XSPEC 的 norm）
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
import struct
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# matplotlib: 强制无 GUI
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# ✅ 你只需要改这里（按你要求）
# =============================================================================
DEFAULT_ROOT = Path("~/data/monk/plot/warmcorona/test/2026.02.04").expanduser().resolve()
DEFAULT_SEARCH_DIR = (DEFAULT_ROOT / "xmm_paramdump_with_error").resolve()
DATA_DIR = (DEFAULT_ROOT / "data").resolve()

# ---- MONK executables ----
SPHERE_EXE = Path("/home/hdw/data/monk/monk_for_rhy/bin/sphere")
SLAB_EXE   = Path("/home/hdw/data/monk/monk_for_rhy/bin/slab")
CALSPEC_EXE = Path("/home/hdw/data/monk/monk_for_rhy/bin/calspec")

# ---- params templates ----
PARAMS_TEMPLATE_SPHERE = Path("~/data/monk/plot/warmcorona/test/test_smooth_10_7/params.txt").expanduser().resolve()
PARAMS_TEMPLATE_SLAB = Path("~/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/params_slab.txt").expanduser().resolve()

# ---- calspec 参数（你要求默认这个）----
CALSPEC_PARAMETER = "-1000.0 0.01 20.0"

# ---- XSPEC energy grid ----
XSPEC_E_MIN = 0.001
XSPEC_E_MAX = 20.0
XSPEC_NBINS = 1000
XSPEC_GRID_MODE = "log"  # "log" or "lin"

# ---- 图的一些默认设置 ----
FLUX_FLOOR = 1e-80


# =============================================================================
# Helpers
# =============================================================================
def die(msg: str):
    raise RuntimeError(msg)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def nonempty_file(p: Path) -> bool:
    return p.is_file() and p.stat().st_size > 0

def read_text_any(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="ignore")

def find_latest_paramdump(search_dir: Path) -> Path:
    if not search_dir.is_dir():
        die(f"[FATAL] Search dir not found: {search_dir}")
    files = sorted(search_dir.glob("xspec_paramdump_with_error_*.txt"),
                   key=lambda x: x.stat().st_mtime, reverse=True)
    if not files:
        die(f"[FATAL] No xspec_paramdump_with_error_*.txt in {search_dir}")
    return files[0]


# =============================================================================
# 1) Parse paramdump -> Tasks
# =============================================================================
HDR_RE = re.compile(
    r"^(?P<src>[^/\s]+)/(?P<obsid>\d{10})\s+rank=(?P<rank>\d+)\s+chi_kind=(?P<kind>\w+)\s+chi_mtime=(?P<mtime>.+)$"
)
SEP_RE = re.compile(r"^=+\s*$")
FLOAT_RE = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"

TE_RE   = re.compile(rf"\bwarmcom\w*te\b.*?\s({FLOAT_RE})\b", re.IGNORECASE)
TAU_RE  = re.compile(rf"\bwarmcom\w*tau\b.*?\s({FLOAT_RE})\b", re.IGNORECASE)
NORM_RE = re.compile(rf"\bwarmcom\w*norm\b.*?\s({FLOAT_RE})\b", re.IGNORECASE)

Z_RE_1 = re.compile(rf"\bz\s*=\s*({FLOAT_RE})\b", re.IGNORECASE)
Z_RE_2 = re.compile(rf"\bwarmcom\w*z\b.*?\s({FLOAT_RE})\b", re.IGNORECASE)

MOD_RE = re.compile(r"atable\{([^}]*warmcom[^}]*\.mod)\}", re.IGNORECASE)

@dataclass
class Task:
    src: str
    obsid: str
    geom: str        # "sphere" or "slab"
    rank: int
    chi_mtime: str
    paramdump: str

    te: Optional[float] = None
    tau: Optional[float] = None
    z: Optional[float] = None
    norm: Optional[float] = None
    warmcom_mod: Optional[str] = None


def split_sections(text: str) -> List[Tuple[Task, List[str]]]:
    lines = text.splitlines()
    out: List[Tuple[Task, List[str]]] = []
    cur_task: Optional[Task] = None
    cur_lines: List[str] = []

    for line in lines:
        m = HDR_RE.match(line.strip())
        if m:
            if cur_task is not None:
                out.append((cur_task, cur_lines))
            cur_task = Task(
                src=m.group("src").strip(),
                obsid=m.group("obsid").strip(),
                geom=m.group("kind").strip().lower(),
                rank=int(m.group("rank")),
                chi_mtime=m.group("mtime").strip(),
                paramdump="",
            )
            cur_lines = []
            continue

        if cur_task is not None and (not SEP_RE.match(line)):
            cur_lines.append(line)

    if cur_task is not None:
        out.append((cur_task, cur_lines))
    return out


def extract_one(section_lines: List[str], regex: re.Pattern) -> Optional[float]:
    for line in section_lines:
        m = regex.search(line)
        if m:
            for g in reversed(m.groups()):
                try:
                    return float(g)
                except Exception:
                    pass
    return None


def extract_mod(section_lines: List[str]) -> Optional[str]:
    for line in section_lines:
        m = MOD_RE.search(line)
        if m:
            return m.group(1).strip()
    return None


def parse_paramdump(paramdump: Path, debug_parse: bool = False) -> List[Task]:
    text = read_text_any(paramdump)
    sections = split_sections(text)

    tasks: List[Task] = []
    for t, sec in sections:
        t.paramdump = str(paramdump)

        t.te = extract_one(sec, TE_RE)
        t.tau = extract_one(sec, TAU_RE)
        t.norm = extract_one(sec, NORM_RE)

        z1 = extract_one(sec, Z_RE_1)
        z2 = extract_one(sec, Z_RE_2)
        t.z = z1 if z1 is not None else z2

        t.warmcom_mod = extract_mod(sec)

        if debug_parse:
            sid = f"{t.src}/{t.obsid} geom={t.geom} rank={t.rank}"
            print(f"\n[DEBUG-parse] {sid}")
            print(f"  te   = {t.te}  | {('OK' if t.te is not None else 'None')}")
            print(f"  tau  = {t.tau}  | {('OK' if t.tau is not None else 'None')}")
            print(f"  z    = {t.z}  | {('OK' if t.z is not None else 'None')}")
            print(f"  norm = {t.norm}  | {('OK' if t.norm is not None else 'None')}")
            print(f"  mod  = {t.warmcom_mod}")

        tasks.append(t)

    return tasks


# =============================================================================
# 2) Directory layout (keep exactly your structure)
# =============================================================================
def task_dir(t: Task) -> Path:
    if t.te is None or t.tau is None:
        te = float("nan"); tau = float("nan")
    else:
        te = float(t.te); tau = float(t.tau)
    return (DATA_DIR / t.src / t.obsid / f"te_{te:.3f}_tau_{tau:.3f}").resolve()

def geom_dir(t: Task) -> Path:
    return task_dir(t) / t.geom

def calspec_dir(t: Task) -> Path:
    return task_dir(t) / "calspec"

def logs_dir(t: Task) -> Path:
    return task_dir(t) / "logs"

def record_json_path(t: Task) -> Path:
    return task_dir(t) / "record.json"

def xspec_qdp_path(t: Task) -> Path:
    return task_dir(t) / "xspec_warmcom_model.qdp"

def xspec_xcm_path(t: Task) -> Path:
    return task_dir(t) / "xspec_warmcom_model.xcm"

def plot_pdf_path(t: Task) -> Path:
    return task_dir(t) / "compare_warmcom_vs_monk.pdf"


# =============================================================================
# 3) MONK runner (skip if already exists)
# =============================================================================
def ensure_env_ok():
    ensure_dir(DATA_DIR)
    if not SPHERE_EXE.is_file():
        die(f"[FATAL] SPHERE_EXE not found: {SPHERE_EXE}")
    if not SLAB_EXE.is_file():
        die(f"[FATAL] SLAB_EXE not found: {SLAB_EXE}")
    if not CALSPEC_EXE.is_file():
        die(f"[FATAL] CALSPEC_EXE not found: {CALSPEC_EXE}")
    if not PARAMS_TEMPLATE_SPHERE.is_file():
        die(f"[FATAL] PARAMS_TEMPLATE_SPHERE not found: {PARAMS_TEMPLATE_SPHERE}")
    if not PARAMS_TEMPLATE_SLAB.is_file():
        die(f"[FATAL] PARAMS_TEMPLATE_SLAB not found: {PARAMS_TEMPLATE_SLAB}")

def monk_outputs_ready(t: Task) -> bool:
    en = calspec_dir(t) / "en.dat"
    fl = calspec_dir(t) / "flux.dat"
    return nonempty_file(en) and nonempty_file(fl)

def write_params_for_geom(t: Task, dst_params: Path):
    if t.te is None or t.tau is None:
        die(f"Missing te/tau for {t.src}_{t.obsid}_{t.geom}")

    if t.geom == "sphere":
        text = read_text_any(PARAMS_TEMPLATE_SPHERE)
        text = re.sub(r"(?m)^\s*te\s*=\s*.*$", f"te = {t.te:.3f}", text, count=1)
        text = re.sub(r"(?m)^\s*tau\s*=\s*.*$", f"tau = {t.tau:.3f}", text, count=1)
        dst_params.write_text(text, encoding="utf-8", errors="ignore")
        return

    if t.geom == "slab":
        text = read_text_any(PARAMS_TEMPLATE_SLAB)
        text = re.sub(r"(?m)^\s*te\s*=\s*.*$", f"te = {t.te:.3f}", text, count=1)
        text = re.sub(r"(?m)^\s*tau\s*=\s*.*$", f"tau = {t.tau:.3f}", text, count=1)
        dst_params.write_text(text, encoding="utf-8", errors="ignore")
        return

    die(f"Unknown geom: {t.geom}")

def run_cmd_show(cmd: str, logfile: Optional[Path] = None, cwd: Optional[Path] = None) -> int:
    """
    运行命令，并把输出显示到终端；如给 logfile，则同时保存（tee）。
    """
    if logfile is not None:
        ensure_dir(logfile.parent)
        full = f"bash -lc \"{cmd} 2>&1 | tee '{logfile}'\""
        p = subprocess.run(full, shell=True, cwd=str(cwd) if cwd else None)
        return p.returncode
    else:
        p = subprocess.run(cmd, shell=True, cwd=str(cwd) if cwd else None)
        return p.returncode

def run_monk_if_needed(t: Task, verbose: bool = False) -> Dict[str, str]:
    ensure_dir(task_dir(t))
    ensure_dir(geom_dir(t))
    ensure_dir(calspec_dir(t))
    ensure_dir(logs_dir(t))

    params_path = geom_dir(t) / "params.txt"
    if not params_path.exists():
        write_params_for_geom(t, params_path)

    info: Dict[str, str] = {}

    if monk_outputs_ready(t):
        if verbose:
            print(f"[SKIP MONK] exists: {t.src}_{t.obsid}_{t.geom} -> {calspec_dir(t)}")
        info["monk_skipped"] = "true"
        return info

    if verbose:
        print(f"[RUN MONK] {t.src}_{t.obsid}_{t.geom} te={t.te} tau={t.tau}")

    exe = SPHERE_EXE if t.geom == "sphere" else SLAB_EXE
    geom_log = logs_dir(t) / f"{t.geom}.log"
    rc = run_cmd_show(f"cd '{geom_dir(t)}' ; sleep 0.01 ; '{exe}'", logfile=geom_log)
    if rc != 0:
        die(f"MONK {t.geom} failed rc={rc} (see {geom_log})")

    cal_log = logs_dir(t) / "calspec.log"
    rel_in = f"../{t.geom}/"
    rc = run_cmd_show(f"cd '{calspec_dir(t)}' && '{CALSPEC_EXE}' {rel_in} {CALSPEC_PARAMETER}", logfile=cal_log)
    if rc != 0:
        die(f"calspec failed rc={rc} (see {cal_log})")

    en = calspec_dir(t) / "en.dat"
    fl = calspec_dir(t) / "flux.dat"
    if (not nonempty_file(en)) or (not nonempty_file(fl)):
        die(f"calspec outputs missing/empty: {en} {fl}")

    info["monk_skipped"] = "false"
    return info


# =============================================================================
# 4) XSPEC export -> QDP (no newpar; sanity check)
# =============================================================================
def parse_showpar_values_from_log(log_text: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    rex_te   = re.compile(r"\bwarmcom\w*te\b\s+(" + FLOAT_RE + r")", re.IGNORECASE)
    rex_tau  = re.compile(r"\bwarmcom\w*tau\b\s+(" + FLOAT_RE + r")", re.IGNORECASE)
    rex_z    = re.compile(r"\bwarmcom\w*z\b\s+(" + FLOAT_RE + r")", re.IGNORECASE)
    rex_norm = re.compile(r"\bwarmcom\w*norm\b\s+(" + FLOAT_RE + r")", re.IGNORECASE)

    for line in log_text.splitlines():
        if "warmcom" not in line.lower():
            continue
        m = rex_te.search(line)
        if m and "te" not in out:
            out["te"] = float(m.group(1))
        m = rex_tau.search(line)
        if m and "tau" not in out:
            out["tau"] = float(m.group(1))
        m = rex_z.search(line)
        if m and "z" not in out:
            out["z"] = float(m.group(1))
        m = rex_norm.search(line)
        if m and "norm" not in out:
            out["norm"] = float(m.group(1))
    return out


def xspec_export_needed(t: Task) -> bool:
    qdp = xspec_qdp_path(t)
    if not nonempty_file(qdp):
        return True

    logp = logs_dir(t) / "xspec_export.log"
    if not nonempty_file(logp):
        return True

    vals = parse_showpar_values_from_log(read_text_any(logp))
    if not all(k in vals for k in ("te", "tau", "z", "norm")):
        return True

    if t.te is None or t.tau is None or t.z is None or t.norm is None:
        return True

    def close(a: float, b: float) -> bool:
        return abs(a - b) <= 1e-6 * max(1.0, abs(b))

    ok = close(vals["te"], float(t.te)) and close(vals["tau"], float(t.tau)) \
         and close(vals["z"], float(t.z)) and close(vals["norm"], float(t.norm))
    return not ok


def run_capture_bytes(cmd_list, input_text: str, cwd: Path, env: Dict[str, str]) -> Tuple[int, str, str]:
    """
    永远以 bytes 捕获 stdout/stderr，自己 decode（replace），避免任何 UnicodeDecodeError。
    返回: (returncode, stdout_text, stderr_text)
    """
    proc = subprocess.run(
        cmd_list,
        input=input_text.encode("utf-8", errors="replace"),
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    stdout_txt = (proc.stdout or b"").decode("utf-8", errors="replace")
    stderr_txt = (proc.stderr or b"").decode("utf-8", errors="replace")
    return proc.returncode, stdout_txt, stderr_txt


def run_xspec_export(t: Task, verbose: bool = False) -> None:
    if t.te is None or t.tau is None or t.z is None or t.norm is None:
        die(f"Missing te/tau/z/norm for XSPEC export: {t.src}_{t.obsid}_{t.geom}")
    if not t.warmcom_mod:
        die(f"Missing warmcom .mod path in paramdump: {t.src}_{t.obsid}_{t.geom}")

    ensure_dir(task_dir(t))
    ensure_dir(logs_dir(t))

    out_qdp = xspec_qdp_path(t)
    out_xcm = xspec_xcm_path(t)
    logp = logs_dir(t) / "xspec_export.log"

    if not xspec_export_needed(t):
        if verbose:
            print(f"[SKIP XSPEC] valid qdp exists: {out_qdp}")
        return

    if out_qdp.exists():
        out_qdp.unlink()

    grid_kw = "log" if XSPEC_GRID_MODE.lower().startswith("log") else "lin"

    cmds = f"""query yes
cpd /null
model atable{{{t.warmcom_mod}}}
{t.te}
{t.tau}
{t.z}
{t.norm}
show par
setplot rebin 1 1
setplot energy
setplot area off
energies {XSPEC_E_MIN} {XSPEC_E_MAX} {XSPEC_NBINS} {grid_kw}
plot model
iplot
wdata {out_qdp.name}
quit
exit
"""

    out_xcm.write_text(cmds, encoding="utf-8", errors="replace")

    env = os.environ.copy()
    env.setdefault("PGPLOT_DEV", "/null")
    env.setdefault("QT_QPA_PLATFORM", "offscreen")

    rc, out_txt, err_txt = run_capture_bytes(["xspec"], cmds, task_dir(t), env)

    logp.write_text(
        "=== STDOUT ===\n" + (out_txt or "") + "\n\n=== STDERR ===\n" + (err_txt or ""),
        encoding="utf-8",
        errors="replace"
    )

    if rc != 0:
        die(f"XSPEC failed rc={rc} (see {logp})")

    if not nonempty_file(out_qdp):
        die(f"XSPEC wdata not produced or empty: {out_qdp} (see {logp})")

    vals = parse_showpar_values_from_log(read_text_any(logp))
    if not all(k in vals for k in ("te", "tau", "z", "norm")):
        out_qdp.unlink(missing_ok=True)
        die(f"XSPEC sanity check failed: cannot parse warmcom params from show par (see {logp})")

    def close(a: float, b: float) -> bool:
        return abs(a - b) <= 1e-6 * max(1.0, abs(b))

    if not (close(vals["te"], float(t.te)) and close(vals["tau"], float(t.tau))
            and close(vals["z"], float(t.z)) and close(vals["norm"], float(t.norm))):
        out_qdp.unlink(missing_ok=True)
        die(f"XSPEC sanity check failed: warmcom params mismatch (see {logp})")

    if verbose:
        print(f"[OK] XSPEC QDP: {out_qdp}")


# =============================================================================
# 5) Read data + Plot
# =============================================================================
def load_dat_file(file_path: str) -> np.ndarray:
    """
    按 little-endian double(<d) 读取二进制 dat 文件。
    """
    with open(file_path, "rb") as f:
        data = f.read()
    if len(data) % 8 != 0:
        raise ValueError(f"{file_path}: 字节长度不是8的整数倍，无法按 double 解析。")
    arr = struct.unpack("<" + "d" * (len(data) // 8), data)
    return np.asarray(arr, dtype=float)

def load_text_1col(file_path: str) -> np.ndarray:
    """
    读取文本格式的一列浮点数（容错：忽略非数字行）。
    """
    xs = []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            try:
                xs.append(float(s.split()[0]))
            except Exception:
                continue
    return np.asarray(xs, dtype=float)

def load_1d_auto(file_path: Path) -> np.ndarray:
    """
    自动判断文件是二进制<double>还是文本：优先尝试二进制，失败再用文本。
    """
    try:
        arr = load_dat_file(str(file_path))
        if arr.size > 0:
            return arr
    except Exception:
        pass

    arr = load_text_1col(str(file_path))
    if arr.size == 0:
        raise ValueError(f"{file_path}: 无法读取（既不像 binary<double> 也不像可解析文本）。")
    return arr

def read_monk_spectrum(t: Task) -> Tuple[np.ndarray, np.ndarray]:
    en_path = calspec_dir(t) / "en.dat"
    fl_path = calspec_dir(t) / "flux.dat"

    en = load_1d_auto(en_path).flatten()
    fl = load_1d_auto(fl_path).flatten()

    n = min(len(en), len(fl))
    en, fl = en[:n], fl[:n]
    m = np.isfinite(en) & np.isfinite(fl) & (en > 0) & (fl > 0)
    return en[m], fl[m]

def read_qdp(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    for line in read_text_any(path).splitlines():
        s = line.strip()
        if (not s) or s[0] in ("!", "@"):
            continue
        parts = s.split()
        try:
            if len(parts) >= 3:
                x = float(parts[0]); y = float(parts[2])
            elif len(parts) >= 2:
                x = float(parts[0]); y = float(parts[1])
            else:
                continue
        except Exception:
            continue
        xs.append(x); ys.append(y)
    if len(xs) < 10:
        die(f"QDP too few points: {path}")
    x = np.array(xs, float)
    y = np.array(ys, float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    return x[m], y[m]

def interp_loglog(x_from: np.ndarray, y_from: np.ndarray, x_to: np.ndarray) -> np.ndarray:
    y = np.full_like(x_to, np.nan, dtype=float)
    m = np.isfinite(x_from) & np.isfinite(y_from) & (x_from > 0) & (y_from > 0)
    if m.sum() < 2:
        return y
    xf = x_from[m]
    yf = y_from[m]
    idx = np.argsort(xf)
    xf = xf[idx]; yf = yf[idx]
    return np.exp(np.interp(np.log(x_to), np.log(xf), np.log(yf), left=np.nan, right=np.nan))

def plot_compare(t: Task, verbose: bool = False) -> None:
    monk_E, monk_F = read_monk_spectrum(t)

    # ✅ 关键：MONK 原始谱乘上 norm 再对比（与你 XSPEC QDP 的 norm 对齐）
    if t.norm is None:
        die(f"Missing norm for scaling MONK spectrum: {t.src}_{t.obsid}_{t.geom}")
    monk_F = monk_F * float(t.norm)

    xspec_E, xspec_F = read_qdp(xspec_qdp_path(t))

    monk_on_x = interp_loglog(monk_E, monk_F, xspec_E)
    ratio = xspec_F / monk_on_x

    fig = plt.figure(figsize=(8.6, 7.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.2], hspace=0.06)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

    ax1.loglog(xspec_E, np.clip(xspec_F, FLUX_FLOOR, np.inf), label="Warmcom (XSPEC atable)")
    ax1.loglog(monk_E, np.clip(monk_F, FLUX_FLOOR, np.inf), label=f"MONK ({t.geom}+calspec) × norm", alpha=0.9)
    ax1.set_ylabel("Flux")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="best", fontsize=9)

    title = f"{t.src}  {t.obsid}  geom={t.geom}\n" \
            f"te={t.te}  tau={t.tau}  z={t.z}  norm={t.norm}"
    ax1.set_title(title, fontsize=10)

    ax2.semilogx(xspec_E, ratio)
    ax2.axhline(1.0, lw=1)
    ax2.set_xlabel("Energy (keV)")
    ax2.set_ylabel("XSPEC / MONK")
    ax2.grid(True, which="both", alpha=0.25)

    rr = ratio[np.isfinite(ratio) & (ratio > 0)]
    if rr.size > 20:
        p1, p99 = np.percentile(rr, [1, 99])
        lo = max(0.2, p1 * 0.8)
        hi = min(5.0, p99 * 1.2)
        if hi > lo:
            ax2.set_ylim(lo, hi)

    out_pdf = plot_pdf_path(t)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    if verbose:
        print(f"[OK] plot: {out_pdf}")


# =============================================================================
# 6) Main
# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paramdump", type=str, default="",
                    help="指定 paramdump；不写则自动选择 DEFAULT_SEARCH_DIR 最新的")
    ap.add_argument("--only-geom", choices=["", "sphere", "slab"], default="",
                    help="只处理某个几何")
    ap.add_argument("--debug-parse", action="store_true", help="打印解析到的 te/tau/z/norm/mod")
    ap.add_argument("--verbose", action="store_true", help="显示更详细日志")
    args = ap.parse_args()

    ensure_env_ok()

    if args.paramdump:
        paramdump = Path(args.paramdump).expanduser().resolve()
        if not paramdump.is_file():
            die(f"[FATAL] --paramdump not found: {paramdump}")
    else:
        paramdump = find_latest_paramdump(DEFAULT_SEARCH_DIR)

    tasks = parse_paramdump(paramdump, debug_parse=args.debug_parse)
    if args.only_geom:
        tasks = [t for t in tasks if t.geom == args.only_geom]

    ensure_dir(DATA_DIR)

    summary = DATA_DIR / "_tasks_summary.txt"
    with summary.open("w", encoding="utf-8") as f:
        f.write(f"[USE] paramdump: {paramdump}\n")
        f.write(f"[USE] data dir : {DATA_DIR}\n")
        f.write(f"[USE] tasks    : {len(tasks)}\n\n")
        for t in tasks:
            f.write(f"{t.geom:6s}  {t.src:18s} {t.obsid}  rank={t.rank}  "
                    f"te={t.te} tau={t.tau} z={t.z} norm={t.norm}\n")

    print(f"[USE] paramdump: {paramdump}")
    print(f"[USE] data dir : {DATA_DIR}")
    print(f"[USE] tasks    : {len(tasks)}")
    print(f"[USE] summary  : {summary}")

    failures = DATA_DIR / "_failures.log"
    ok, fail = 0, 0
    t0 = time.time()

    for t in tasks:
        tag = f"{t.src}_{t.obsid}_{t.geom}"
        try:
            if t.te is None or t.tau is None or t.z is None or t.norm is None:
                die(f"Missing te/tau/z/norm in paramdump for task {tag}")
            if not t.warmcom_mod:
                die(f"Missing warmcom .mod path in paramdump for task {tag}")

            monk_info = run_monk_if_needed(t, verbose=args.verbose)
            run_xspec_export(t, verbose=args.verbose)
            plot_compare(t, verbose=args.verbose)

            rec = asdict(t)
            rec.update({
                "data_dir": str(task_dir(t)),
                "geom_dir": str(geom_dir(t)),
                "calspec_dir": str(calspec_dir(t)),
                "logs_dir": str(logs_dir(t)),
                "monk_skipped": monk_info.get("monk_skipped", "unknown"),
                "calspec_parameter": CALSPEC_PARAMETER,
                "xspec_qdp": str(xspec_qdp_path(t)),
                "xspec_xcm": str(xspec_xcm_path(t)),
                "plot_pdf": str(plot_pdf_path(t)),
                "monk_flux_scaled_by_norm": True,
            })
            record_json_path(t).write_text(json.dumps(rec, ensure_ascii=False, indent=2),
                                           encoding="utf-8", errors="replace")

            ok += 1

        except Exception as e:
            fail += 1
            with failures.open("a", encoding="utf-8") as ff:
                ff.write(f"\n[{tag}] ERROR:\n{e}\n")
            print(f"[FAIL] {tag}: {e}", flush=True)

    dt = time.time() - t0
    print(f"\n[DONE] ok={ok} fail={fail} wall={dt:.1f}s")
    if fail:
        print(f"[DONE] failures: {failures}")


if __name__ == "__main__":
    main()
