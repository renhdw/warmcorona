#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
warmcom_validate_3models_heatmaps.py

✅ err = |1 - MODEL/REF|
✅ 无效点( MODEL<=0 or REF<=0 or non-finite ) -> NaN（不画）
✅ 尾巴处理：
   若 MODEL 在某处首次 <= max(MODEL_this_row)*RATIO_TAIL，则从该处起：
     - 若该点本来是“有效误差”(finite)，则 err = 0（仍画出来）
     - 若该点本来无效(NaN)，保持 NaN（不画）
✅ 每个模型只启动一次 XSPEC（共 3 次），XSPEC 内部循环 N 点
✅ wdata 放到 iplot：iplot -> wdata -> quit；最后才 exit XSPEC

✅ 统一色阶（固定，不随每次运行改变）：
   - 0–0.05：渐变（用于分辨小误差）
   - >=0.05：全部同一个“深蓝色”（0.05–50 与 >50 不再区分）
   - NaN：透明/白（不画）

✅ REF 定义（warmcom 版本）：
   - 第一张图（FULL=warmcomsphere）：REF = FULL（err = |1 - FULL/FULL| ~ 0，仅 sanity check）
   - 其余两张图（tetest/tautest）：REF = FULL（err = |1 - MODEL/FULL|）

✅ 每页下方加一张散点图（不增加页数）
   - y：Selected point No
   - x：该行误差（compute_error 之后）最大值 max(err)

✅ 新增：每个模型额外再输出 1 页 “All-points ratio-vs-error 散点图（rasterized）”
   - 横坐标：F_MODEL(E) / max_E F_MODEL(E)
   - 纵坐标：|1 - MODEL/REF|（compute_error 后）
   - 所有点合并，不区分组号/No
   - scatter rasterized=True 防止 PDF 体积爆
   - 点数过大自动下采样（--max-points）
"""

import os
import re
import time
import shutil
import argparse
import subprocess
import numpy as np
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, ListedColormap
from matplotlib.ticker import LogLocator, LogFormatterMathtext, FuncFormatter

# =========================
# Defaults (YOU)
# =========================
OUT_DIR_DEFAULT = os.path.expanduser("~/data/monk/plot/warmcorona/test/2026.01.26/warmcom_validate_3models_heatmaps_sphere_v2/log")

MODEL_BASEDIR = os.path.expanduser("~/data/monk/plot/warmcorona/model/warmcom/warmcom_sphere_log")
MODEL_DIR_FULL      = os.path.join(MODEL_BASEDIR, "warmcom_sphere")
MODEL_DIR_TETEST    = os.path.join(MODEL_BASEDIR, "warmcom_sphere_tetest")
MODEL_DIR_TAUTEST   = os.path.join(MODEL_BASEDIR, "warmcom_sphere_tautest")

PKG_NAME_FULL     = "warmcom_sphere"
MODEL_NAME_FULL   = "warmcomsphere"

PKG_NAME_TETEST   = "warmcom_sphere_tetest"
MODEL_NAME_TETEST = "warmcomspheretetest"

PKG_NAME_TAUTEST  = "warmcom_sphere_tautest"
MODEL_NAME_TAUTEST= "warmcomspheretautest"

# warmcom base grid dir (where chose.log and te_*/tau_* dirs live)
WARMCOM_BASEDIR_DEFAULT = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7/_exports"
WARMCOM_BASEDIR_DEFAULT = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7/_exports"
CHOSE_LOG_SUBPATH = "chose.log"

# XSPEC plotting energy grid
E_MIN = 0.001
E_MAX = 20.0
NBINS = 1000
E_GRID_MODE = "log"

# Model params defaults
Z_DEFAULT = 0.0
NORM_DEFAULT = 1.0

# Fixed colormap config
FIXED_VMAX = 50.0
THRESH_SOLID = 0.05
RATIO_TAIL = 1e-20
SOLID_BLUE = "#08306b"

# All-points scatter: x_min keep consistent with tail ratio
X_RATIO_MIN = RATIO_TAIL

# All-points scatter: for log-y display of err=0
Y_ERR_FLOOR = 1e-6

# =========================
# Helpers
# =========================
def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def read_qdp_model(path: str):
    xs, ys = [], []
    with open(path, "r") as f:
        for line in f:
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
    if not xs:
        raise RuntimeError(f"QDP 为空：{path}")
    return np.array(xs, float), np.array(ys, float)

def _tail_zero_by_model_threshold_keep_nan(err: np.ndarray, y_model: np.ndarray, ratio: float = RATIO_TAIL) -> np.ndarray:
    err = np.asarray(err, float).copy()
    y_model = np.asarray(y_model, float)

    m_ok = np.isfinite(y_model) & (y_model > 0)
    if not np.any(m_ok):
        return err

    m_max = float(np.nanmax(y_model[m_ok]))
    if (not np.isfinite(m_max)) or (m_max <= 0.0):
        return err

    thr = m_max * float(ratio)
    idx_candidates = np.flatnonzero(m_ok & (y_model <= thr))
    if idx_candidates.size == 0:
        return err

    cut = int(idx_candidates[0])
    tail = np.arange(err.size) >= cut
    err[tail & np.isfinite(err)] = 0.0
    return err

def compute_error(y_model: np.ndarray, y_ref: np.ndarray) -> np.ndarray:
    y_model = np.asarray(y_model, float)
    y_ref   = np.asarray(y_ref, float)

    ok = np.isfinite(y_model) & np.isfinite(y_ref) & (y_model > 0) & (y_ref > 0)
    err = np.full_like(y_model, np.nan, dtype=float)
    err[ok] = np.abs(1.0 - (y_model[ok] / y_ref[ok]))
    err = _tail_zero_by_model_threshold_keep_nan(err, y_model, ratio=RATIO_TAIL)
    return err

def make_edges_from_centers_log(centers: np.ndarray) -> np.ndarray:
    c = np.asarray(centers, float)
    if np.any(c <= 0):
        raise ValueError("Energy centers must be >0 for log edges.")
    edges = np.empty(len(c) + 1, dtype=float)
    mid = np.sqrt(c[:-1] * c[1:])
    edges[1:-1] = mid
    edges[0]  = c[0] * (c[0] / mid[0])
    edges[-1] = c[-1] * (c[-1] / mid[-1])
    return edges

def make_levels_solid_after_threshold(vmax: float = FIXED_VMAX, thr: float = THRESH_SOLID):
    thr = float(thr); vmax = float(vmax)
    if not (thr > 0 and vmax > thr):
        raise ValueError("Need vmax>thr>0")

    fine = np.array([
        0.0,
        1e-6, 2e-6, 5e-6,
        1e-5, 2e-5, 5e-5,
        1e-4, 2e-4, 5e-4,
        1e-3, 2e-3, 5e-3,
        1e-2, 2e-2, 3e-2, 4e-2,
        thr
    ], dtype=float)

    fine = fine[(fine >= 0) & (fine <= thr)]
    if fine[-1] != thr:
        fine = np.r_[fine, thr]

    levels = np.r_[fine, vmax]
    if not np.all(np.diff(levels) > 0):
        raise RuntimeError("levels must be strictly increasing.")
    return levels

def make_cmap_solid_after_threshold(levels: np.ndarray, thr: float = THRESH_SOLID, solid_color: str = SOLID_BLUE):
    levels = np.asarray(levels, float)
    nbins = len(levels) - 1
    if nbins <= 1:
        raise ValueError("levels too short")

    thr = float(thr)
    j = np.where(np.isclose(levels, thr, rtol=0.0, atol=0.0))[0]
    if j.size == 0:
        raise RuntimeError("thr must be exactly included in levels.")
    j = int(j[0])
    n_grad = j
    n_solid = nbins - n_grad
    if n_grad <= 0 or n_solid <= 0:
        raise RuntimeError("Need both gradient bins and solid bins.")

    cmap_grad = LinearSegmentedColormap.from_list(
        "BeigeToLightRed_0toThr",
        ["#f2e6c9", "#f6e1cf", "#f7d6c8", "#f6c7bd", "#f4b8a5"],
        N=512
    )
    colors_grad = cmap_grad(np.linspace(0, 1, n_grad, endpoint=True))

    rgba_solid = matplotlib.colors.to_rgba(solid_color)
    colors_solid = np.tile(np.array(rgba_solid)[None, :], (n_solid, 1))

    colors = np.vstack([colors_grad, colors_solid])
    cmap = ListedColormap(colors, name="GradThenSolidBlue")
    cmap.set_bad(color="white", alpha=0.0)  # NaN blank/transparent
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    return norm, cmap

def row_max_error_value(err_mat: np.ndarray) -> np.ndarray:
    err_mat = np.asarray(err_mat, float)
    with np.errstate(all="ignore"):
        mx = np.nanmax(err_mat, axis=1)
    return mx

def plot_error_heatmap_with_rowmax(pdf: PdfPages, E: np.ndarray, err_mat: np.ndarray,
                                   y_labels: list, title: str, subtitle: str,
                                   norm, cmap, thr: float = THRESH_SOLID):
    Ny, _ = err_mat.shape
    E_edges = make_edges_from_centers_log(E)
    y_edges = np.arange(Ny + 1)
    y_center = np.arange(Ny) + 0.5

    err_plot = np.ma.masked_invalid(err_mat)

    row_max = row_max_error_value(err_mat)
    ok = np.isfinite(row_max)

    fig = plt.figure(figsize=(11.5, 8.2))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4.6, 1.6], hspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    pcm = ax.pcolormesh(E_edges, y_edges, err_plot, shading="auto",
                        cmap=cmap, norm=norm, rasterized=True)

    ax.set_xscale("log")
    ax.set_xlim(E_edges[0], E_edges[-1])
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
    ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))

    ax.set_xlabel("Energy [keV]")
    ax.set_ylabel("Selected point No")

    if Ny <= 20:
        ax.set_yticks(np.arange(Ny) + 0.5)
        ax.set_yticklabels(y_labels)
    else:
        yt = np.linspace(0, Ny - 1, 10).astype(int)
        ax.set_yticks(yt + 0.5)
        ax.set_yticklabels([y_labels[i] for i in yt])

    ax.set_title(title + "\n" + subtitle, fontsize=11)

    def _cbar_fmt(x, pos):
        if x == 0:
            return "0"
        if x < 1e-3:
            return f"{x:.0e}"
        if x < 0.1:
            return f"{x:.3f}"
        if x < 1:
            return f"{x:.2f}"
        if x < 10:
            return f"{x:.1f}"
        return f"{x:.0f}"

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label(
        rf"Error  $|1-\mathrm{{MODEL}}/\mathrm{{REF}}|$  (0–{thr:g} gradient; >= {thr:g} solid blue; NaN blank)",
        fontsize=10
    )
    ticks = [0, 1e-4, 1e-3, 1e-2, 0.02, 0.03, 0.04, thr, 0.1, 1, 10, FIXED_VMAX]
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(_cbar_fmt))

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.scatter(row_max[ok], y_center[ok], s=18, rasterized=True)
    ax2.set_xscale("log")
    ax2.set_xlim(1e-3, FIXED_VMAX)
    ax2.set_ylim(0, Ny)
    ax2.set_ylabel("Selected point No")
    ax2.set_xlabel(r"Max error per row  $\max\left(|1-\mathrm{MODEL}/\mathrm{REF}|\right)$")

    if Ny <= 20:
        ax2.set_yticks(np.arange(Ny) + 0.5)
        ax2.set_yticklabels(y_labels)
    else:
        yt = np.linspace(0, Ny - 1, 10).astype(int)
        ax2.set_yticks(yt + 0.5)
        ax2.set_yticklabels([y_labels[i] for i in yt])

    ax2.grid(True, which="both", axis="x", alpha=0.3)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

# =========================
# NEW: per-model ALLPOINTS ratio-vs-error (rasterized scatter)
# =========================
def plot_allpoints_ratio_vs_error(pdf: PdfPages,
                                  flux_mat: np.ndarray,
                                  err_mat: np.ndarray,
                                  model_label: str,
                                  subtitle: str,
                                  x_min: float = X_RATIO_MIN,
                                  y_floor: float = Y_ERR_FLOOR,
                                  max_points: int = 4_000_000,
                                  seed: int = 123,
                                  s: float = 1.0,
                                  alpha: float = 0.15):
    """
    一页：所有点合并的 ratio-vs-error 散点图（rasterized）
      x = F(E)/max_E(F)
      y = err(E)

    注意：
      - err 使用 compute_error 输出（含 NaN / tail->0）
      - y=0 为了 log-y 显示，用 y_floor 替换（仅画图）
      - 点太多自动下采样（max_points）
    """
    F = np.asarray(flux_mat, float)
    E = np.asarray(err_mat, float)
    if F.shape != E.shape:
        raise RuntimeError("flux_mat and err_mat shape mismatch.")

    F_pos = np.where(np.isfinite(F) & (F > 0), F, np.nan)
    with np.errstate(all="ignore"):
        row_max = np.nanmax(F_pos, axis=1)
    row_max[~np.isfinite(row_max)] = np.nan

    ratio_mat = F / row_max[:, None]

    ok = np.isfinite(E) & np.isfinite(ratio_mat) & (ratio_mat > 0)
    if not np.any(ok):
        fig = plt.figure(figsize=(9.6, 6.6))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No valid points for ratio-error scatter.",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title(f"All-points ratio vs error: {model_label}\n{subtitle}", fontsize=11)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        return

    x = ratio_mat[ok]
    y = E[ok]

    x = np.clip(x, x_min, 1.0)
    y = np.clip(y, y_floor, None)

    n_raw = x.size
    if (max_points is not None) and (n_raw > int(max_points)):
        rng = np.random.default_rng(seed)
        idx = rng.choice(n_raw, size=int(max_points), replace=False)
        x = x[idx]
        y = y[idx]
    n_plot = x.size

    fig = plt.figure(figsize=(9.6, 6.6))
    ax = fig.add_subplot(111)

    ax.set_rasterization_zorder(0)
    ax.scatter(x, y, s=s, alpha=alpha, linewidths=0, rasterized=True, zorder=0)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_min, 1.05)
    ax.set_ylim(y_floor, FIXED_VMAX)

    ax.grid(True, which="both", alpha=0.25)
    ax.set_xlabel(r"$F_{\rm MODEL}(E)\ /\ \max_E F_{\rm MODEL}(E)$")
    ax.set_ylabel(r"$|1-\mathrm{MODEL}/\mathrm{REF}|$  (0 shown as floor for log-y)")

    ax.set_title(
        f"All-points ratio vs error (rasterized): {model_label}\n{subtitle}\n"
        f"N(raw)={n_raw:,} ; N(plotted)={n_plot:,} ; x_min={x_min:g}(=RATIO_TAIL) ; y_floor={y_floor:g}",
        fontsize=11
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

# =========================
# Parse chose.log to get available (te,tau)
# =========================
def parse_chose_log(base_dir: str):
    logpath = os.path.join(base_dir, CHOSE_LOG_SUBPATH)
    if not os.path.exists(logpath):
        raise RuntimeError(f"Cannot find chose.log: {logpath}")

    pat = re.compile(r"te_([0-9]*\.?[0-9]+)_tau_([0-9]*\.?[0-9]+)")
    pairs = []
    with open(logpath, "r", encoding="utf-8", errors="ignore") as f:
        for ln, line in enumerate(f, start=1):
            if ln <= 3:
                continue
            m = pat.search(line)
            if not m:
                continue
            te = float(m.group(1))
            tau = float(m.group(2))
            pairs.append((te, tau))

    if not pairs:
        raise RuntimeError("No (te,tau) parsed from chose.log")

    pairs = sorted(list(set(pairs)), key=lambda x: (x[0], x[1]))
    return pairs

# =========================
# XSPEC batch per model (ONE session, loop N jobs)
# warmcom params: te, tau, z, norm
# =========================
def xspec_batch_write_qdp_one_model_warmcom(
    model_dir: str, pkg: str, model_name: str,
    jobs: list,
    workdir: str,
    out_qdp_dir: str,
    emin: float, emax: float, nbins: int, grid_mode: str,
    out_xcm: str,
    warmcom_basedir: str
):
    if not jobs:
        raise ValueError("jobs empty")

    ensure_dir(workdir)
    ensure_dir(out_qdp_dir)
    ensure_dir(os.path.dirname(out_xcm))

    grid_kw = "log" if str(grid_mode).lower().startswith("log") else "lin"

    p0 = jobs[0]["params"]
    ptxt0 = "\n".join(str(x) for x in p0)

    for job in jobs:
        No = int(job["No"])
        qdp_base = f"{model_name}_No{No}.qdp"
        for fp in (os.path.join(workdir, qdp_base), os.path.join(out_qdp_dir, qdp_base)):
            if os.path.exists(fp):
                try:
                    os.remove(fp)
                except Exception:
                    pass

    lines = []
    lines.append("query yes")
    lines.append("cpd /null")
    lines.append(f"lmod {pkg} {model_dir}")
    lines.append(f"model {model_name}")
    lines.append(ptxt0)

    lines.append("setplot rebin 1 1")
    lines.append("setplot energy")
    lines.append("setplot area off")
    lines.append(f"energies {emin} {emax} {nbins} {grid_kw}")

    for job in jobs:
        No = int(job["No"])
        te, tau, z, norm = job["params"]
        qdp_base = f"{model_name}_No{No}.qdp"

        lines.append(f"# --- No={No} ---")
        lines.append(f"newpar 1 {te}")
        lines.append(f"newpar 2 {tau}")
        lines.append(f"newpar 3 {z}")
        lines.append(f"newpar 4 {norm}")
        lines.append("plot model")

        lines.append("iplot")
        lines.append(f"wdata {qdp_base}")
        lines.append("quit")

    lines.append("exit")

    cmds = "\n".join(lines) + "\n"
    with open(out_xcm, "w", encoding="utf-8") as f:
        f.write(cmds)

    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env.setdefault("PGPLOT_DEV", "/null")
    env["WARMCOM_BASEDIR"] = warmcom_basedir

    proc = subprocess.run(
        ["xspec"],
        input=cmds,
        text=True,
        cwd=workdir,
        env=env,
        capture_output=True
    )
    if proc.returncode != 0:
        print("\n[XSPEC STDOUT tail]\n", proc.stdout[-5000:])
        print("\n[XSPEC STDERR tail]\n", proc.stderr[-5000:])
        raise RuntimeError(f"XSPEC batch(QDP) failed for model={model_name}: returncode={proc.returncode}")

    missing = []
    for job in jobs:
        No = int(job["No"])
        qdp_base = f"{model_name}_No{No}.qdp"
        qdp_abs  = os.path.join(workdir, qdp_base)
        qdp_dst  = os.path.join(out_qdp_dir, qdp_base)
        if (not os.path.exists(qdp_abs)) or os.path.getsize(qdp_abs) == 0:
            missing.append(qdp_base)
            continue
        shutil.move(qdp_abs, qdp_dst)

    if missing:
        print("\n[XSPEC STDOUT tail]\n", proc.stdout[-5000:])
        print("\n[XSPEC STDERR tail]\n", proc.stderr[-5000:])
        raise RuntimeError(
            f"Missing/empty QDP files for model={model_name}: "
            f"{missing[:8]}{'...' if len(missing) > 8 else ''}"
        )

# =========================
# Main
# =========================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default=OUT_DIR_DEFAULT, help="输出目录")
    parser.add_argument("-n", "--ntrial", type=int, default=10, help="随机选点数量（默认10）")
    parser.add_argument("--seed", type=int, default=None, help="随机种子（默认None，每次随机）")

    parser.add_argument("--basedir", type=str, default=WARMCOM_BASEDIR_DEFAULT,
                        help="WARMCOM_BASEDIR (contains chose.log and te_*/tau_* dirs)")
    parser.add_argument("--z", type=float, default=Z_DEFAULT, help="redshift z")
    parser.add_argument("--norm", type=float, default=NORM_DEFAULT, help="XSPEC norm (model ignores it)")

    # allpoints scatter control
    parser.add_argument("--max-points", type=int, default=4_000_000,
                        help="All-points scatter: 最大点数（超过则随机下采样）")
    parser.add_argument("--scatter-alpha", type=float, default=0.15,
                        help="All-points scatter: alpha")
    parser.add_argument("--scatter-size", type=float, default=1.0,
                        help="All-points scatter: marker size")

    parser.add_argument("--keep-temp", action="store_true", help="保留中间 qdp/xcm（默认删除）")
    args = parser.parse_args()

    t0 = time.time()
    ensure_dir(args.outdir)

    run_log_txt   = os.path.join(args.outdir, "run_log.txt")
    out_pdf       = os.path.join(args.outdir, "error_heatmaps_3models_with_allpoints_scatter.pdf")
    selected_txt  = os.path.join(args.outdir, "selected_points_te_tau.txt")

    tmp_dir  = os.path.join(args.outdir, "_tmp_xspec_batch3")
    work_dir = os.path.join(tmp_dir, "_work")
    qdp_dir  = os.path.join(tmp_dir, "_qdp_outputs")
    xcm_dir  = os.path.join(tmp_dir, "_xcm_scripts")
    ensure_dir(work_dir); ensure_dir(qdp_dir); ensure_dir(xcm_dir)

    print("== Parse chose.log to get available (te,tau) ==")
    pairs = parse_chose_log(args.basedir)
    print(f"  available points = {len(pairs)}")
    if args.ntrial <= 0 or args.ntrial > len(pairs):
        raise RuntimeError(f"ntrial must be in [1, {len(pairs)}]")

    rng = np.random.default_rng(args.seed)
    idx = np.sort(rng.choice(len(pairs), size=args.ntrial, replace=False))
    sel_pairs = [pairs[i] for i in idx]

    with open(selected_txt, "w", encoding="utf-8") as f:
        f.write("# No  te  tau\n")
        for ii, (te, tau) in enumerate(sel_pairs, start=1):
            f.write(f"{ii:4d}  {te:.6g}  {tau:.6g}\n")
    print("  ->", selected_txt)

    jobs = []
    for ii, (te, tau) in enumerate(sel_pairs, start=1):
        jobs.append({"No": ii, "params": [float(te), float(tau), float(args.z), float(args.norm)]})

    models = [
        ("warmcomsphere(FULL)",   MODEL_DIR_FULL,   PKG_NAME_FULL,   MODEL_NAME_FULL),
        ("warmcomspheretetest",   MODEL_DIR_TETEST, PKG_NAME_TETEST, MODEL_NAME_TETEST),
        ("warmcomspheretautest",  MODEL_DIR_TAUTEST,PKG_NAME_TAUTEST,MODEL_NAME_TAUTEST),
    ]

    print("\n== Run XSPEC (QDP batch): 3 models => 3 XSPEC sessions ==")
    for model_label, model_dir, pkg, model_name in models:
        print(f"\n--- Model: {model_label} ---")
        out_xcm = os.path.join(xcm_dir, f"batch_{model_name}_{len(jobs)}pts.xcm")
        xspec_batch_write_qdp_one_model_warmcom(
            model_dir=model_dir, pkg=pkg, model_name=model_name,
            jobs=jobs,
            workdir=work_dir,
            out_qdp_dir=qdp_dir,
            emin=E_MIN, emax=E_MAX, nbins=NBINS, grid_mode=E_GRID_MODE,
            out_xcm=out_xcm,
            warmcom_basedir=args.basedir
        )
        print(f"  QDP collected in: {qdp_dir}")

    # -------------------------
    # Parse FULL first
    # -------------------------
    print("\n== Parse FULL QDP first (build FULL as REF) ==")
    E_ref = None
    Ffull_mat = []
    for job in jobs:
        No = int(job["No"])
        qdp_path = os.path.join(qdp_dir, f"{MODEL_NAME_FULL}_No{No}.qdp")
        Ex, Fm = read_qdp_model(qdp_path)
        if E_ref is None:
            E_ref = Ex.copy()
        else:
            if len(Ex) != len(E_ref) or (not np.allclose(Ex, E_ref, rtol=0.0, atol=0.0)):
                raise RuntimeError(f"Energy grid mismatch in FULL No={No}")
        Ffull_mat.append(Fm)
    Ffull_mat = np.vstack(Ffull_mat)

    # compute error mats + collect flux mats for all models
    err_mats = {}
    flux_mats = {}
    flux_mats["warmcomsphere(FULL)"] = Ffull_mat

    print("\n== Compute error for FULL: err = |1 - FULL/FULL| ==")
    err_full = []
    for i in range(len(jobs)):
        err_full.append(compute_error(Ffull_mat[i], Ffull_mat[i]))
    err_mats["warmcomsphere(FULL)"] = np.vstack(err_full)

    for model_label, _mdir, _pkg, model_name in models[1:]:
        print(f"\n== Parse QDP + compute errors: {model_label}  (err = |1 - MODEL/FULL|) ==")
        F_rows = []
        err_rows = []
        for i, job in enumerate(jobs):
            No = int(job["No"])
            qdp_path = os.path.join(qdp_dir, f"{model_name}_No{No}.qdp")
            Ex, Fm = read_qdp_model(qdp_path)
            if len(Ex) != len(E_ref) or (not np.allclose(Ex, E_ref, rtol=0.0, atol=0.0)):
                raise RuntimeError(f"Energy grid mismatch: {model_name} No={No}")
            F_rows.append(Fm)
            err_rows.append(compute_error(Fm, Ffull_mat[i]))
        flux_mats[model_label] = np.vstack(F_rows)
        err_mats[model_label] = np.vstack(err_rows)

    # colormap/norm
    levels = make_levels_solid_after_threshold(vmax=FIXED_VMAX, thr=THRESH_SOLID)
    fixed_norm, fixed_cmap = make_cmap_solid_after_threshold(levels, thr=THRESH_SOLID, solid_color=SOLID_BLUE)

    subtitle = (
        f"Selected N={args.ntrial}, seed={args.seed}; basedir={args.basedir}\n"
        f"z={args.z:g}, norm={args.norm:g} (note: model ignores norm)\n"
        f"Error = |1 - MODEL/REF| ; invalid: MODEL<=0 or REF<=0 or non-finite -> NaN(blank) ; "
        f"tail: MODEL<=max(MODEL)*{RATIO_TAIL:g} -> err=0 if finite (NaN stays NaN)\n"
        f"REF rule: FULL uses FULL; others use FULL ; "
        f"colormap: 0–{THRESH_SOLID:g} gradient, >= {THRESH_SOLID:g} solid blue\n"
        f"Allpoints scatter: x=F/max(F), x_min={X_RATIO_MIN:g}(=RATIO_TAIL), y_floor={Y_ERR_FLOOR:g}"
    )

    y_labels = [str(int(job["No"])) for job in jobs]

    print("\n== Write PDF: 3 heatmaps + 3 allpoints scatters ==")
    with PdfPages(out_pdf) as pdf:
        # 3 pages heatmaps (each with bottom row-max scatter)
        for model_label, _, _, _ in models:
            plot_error_heatmap_with_rowmax(
                pdf,
                E=E_ref,
                err_mat=err_mats[model_label],
                y_labels=y_labels,
                title=f"Error Heatmap: {model_label}",
                subtitle=subtitle,
                norm=fixed_norm,
                cmap=fixed_cmap,
                thr=THRESH_SOLID
            )

        # + 3 pages allpoints ratio-vs-error rasterized scatter
        for model_label, _, _, _ in models:
            plot_allpoints_ratio_vs_error(
                pdf,
                flux_mat=flux_mats[model_label],
                err_mat=err_mats[model_label],
                model_label=model_label,
                subtitle=subtitle,
                x_min=X_RATIO_MIN,
                y_floor=Y_ERR_FLOOR,
                max_points=args.max_points,
                seed=(args.seed if args.seed is not None else 123),
                s=args.scatter_size,
                alpha=args.scatter_alpha
            )

    print("  ->", out_pdf)

    if not args.keep_temp:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        print("\n== Temp files removed ==")
    else:
        print("\n== Temp files kept ==")
        print("  tmp_dir =", tmp_dir)

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    dt = time.time() - t0
    with open(run_log_txt, "w", encoding="utf-8") as f:
        f.write(f"timestamp: {ts}\n")
        f.write(f"outdir: {args.outdir}\n")
        f.write(f"basedir: {args.basedir}\n")
        f.write(f"ntrial: {args.ntrial}\n")
        f.write(f"seed: {args.seed}\n")
        f.write(f"z: {args.z}\n")
        f.write(f"norm: {args.norm}\n")
        f.write(f"pdf: {out_pdf}\n")
        f.write(f"selected_points: {selected_txt}\n")
        f.write(f"keep_temp: {args.keep_temp}\n")
        f.write(f"total_time_sec: {dt:.3f}\n")
        f.write(f"solid_threshold: {THRESH_SOLID}\n")
        f.write(f"solid_color: {SOLID_BLUE}\n")
        f.write(f"fixed_vmax: {FIXED_VMAX}\n")
        f.write(f"tail_ratio: {RATIO_TAIL}\n")
        f.write(f"x_ratio_min: {X_RATIO_MIN} (must equal tail_ratio)\n")
        f.write(f"y_err_floor: {Y_ERR_FLOOR}\n")
        f.write(f"allpoints_max_points: {args.max_points}\n")
        f.write(f"allpoints_scatter_alpha: {args.scatter_alpha}\n")
        f.write(f"allpoints_scatter_size: {args.scatter_size}\n")
        f.write("invalid_rule: MODEL<=0 or REF<=0 or non-finite -> NaN(blank)\n")
        f.write("tail_rule: if MODEL <= max(MODEL_this_row)*tail_ratio then err=0 for finite err; NaN stays NaN\n")
        f.write("ref_rule: FULL uses FULL; others use FULL\n")
        f.write("cmap: 0–thr gradient; >=thr solid blue; NaN transparent\n")
        f.write("levels:\n")
        for x in levels:
            f.write(f"{x:.16g}\n")

    print("  ->", run_log_txt)
    print("\n✅ DONE")
    print(f"  total time = {dt:.2f} s")
    print("\nRun examples:")
    print("  python3 warmcom_validate_3models_heatmaps.py -n 100")
    print("  python3 warmcom_validate_3models_heatmaps.py -n 2000 --max-points 800000")
    print("  python3 warmcom_validate_3models_heatmaps.py -n 2000 --max-points 200000 --scatter-alpha 0.08")

if __name__ == "__main__":
    main()
