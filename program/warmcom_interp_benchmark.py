#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
warmcom_interp_benchmark.py

✅ 读取 warmcom 导出的二进制 .dat (little-endian double)
✅ 横坐标严格使用 en.dat（来自每个点目录）
✅ 固定一个维度（--mode te 或 --mode tau），在另一个维度上做 1D 插值误差测试：
   - slice 上所有点按参数升序
   - 偶数索引为网格(solid)，奇数索引为测试(dashed)
✅ 四种插值：
   1) linp_linF
   2) linp_logF
   3) logp_linF
   4) logp_logF

✅ 兼容两种目录布局（自动探测）：
  A) BASEDIR/te_*/tau_*/calspec/en.dat   (默认假设)
  B) BASEDIR/calspec/te_*/tau_*/en.dat

输出：
- PDF：1页谱线示意 + 4页误差热图
- summary.txt：统计
- used_dirs.txt：本次实际使用到的真实目录（可复现/排错）
"""

import os
import re
import struct
import argparse
import numpy as np
from typing import Dict, List, Tuple, Optional
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# =====================
# DEFAULTS
# =====================
BASEDIR_DEFAULT = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"

OUTBASE_DEFAULT = os.path.expanduser("~/data/monk/plot/warmcorona/test/2026.01.26")
OUTSUB_DEFAULT  = "warmcom_interp_benchmark"

# extra layer
POINT_SUBDIR_DEFAULT = "calspec"

# default filenames inside each point leaf dir
EN_FILE_DEFAULT   = "en.dat"
FLUX_FILE_DEFAULT = "flux.dat"

# slice match tolerance
ATOL_DEFAULT = 1e-10

# numeric safety
FLOOR_POS = 1e-300

# plot styles
LINE_ALPHA_SOLID = 0.55
LINE_ALPHA_DASH  = 0.75
LINE_W_SOLID = 0.9
LINE_W_DASH  = 1.0

# heatmap clip
ERR_CLIP = 5.0


# =====================
# utils
# =====================
def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def load_dat_file(file_path: str) -> np.ndarray:
    """读取二进制 double array（little-endian）"""
    with open(file_path, "rb") as f:
        data = f.read()
    if len(data) % 8 != 0:
        raise ValueError(f"{file_path}: 字节长度不是8的整数倍，无法按 double 解析。")
    n = len(data) // 8
    arr = struct.unpack("<" + "d" * n, data)
    return np.asarray(arr, dtype=float)

def scan_te_tau_dirs(root: str) -> List[Tuple[float, float, str]]:
    """
    扫描 root 下真实存在的 te_*_tau_* 子目录，返回：
      [(te, tau, dirname), ...]
    dirname 为原样目录名（避免 0.10 vs 0.1）
    """
    pat = re.compile(r"^te_([0-9]*\.?[0-9]+)_tau_([0-9]*\.?[0-9]+)$")
    out: List[Tuple[float, float, str]] = []
    if not os.path.isdir(root):
        return out

    for name in os.listdir(root):
        p = os.path.join(root, name)
        if not os.path.isdir(p):
            continue
        m = pat.match(name)
        if not m:
            continue
        te = float(m.group(1))
        tau = float(m.group(2))
        out.append((te, tau, name))

    if not out:
        return out

    # 去重（同一个 (te,tau) 多目录，保留第一个）
    uniq: Dict[Tuple[float, float], str] = {}
    for te, tau, name in out:
        key = (te, tau)
        if key not in uniq:
            uniq[key] = name

    out2 = [(te, tau, uniq[(te, tau)]) for (te, tau) in uniq.keys()]
    out2.sort(key=lambda x: (x[0], x[1]))
    return out2

def detect_layout(basedir: str, point_subdir: str,
                  en_file: str, flux_file: str) -> Tuple[str, List[Tuple[float, float, str]]]:
    """
    自动探测布局：
      A) basedir/te_*/tau_*/point_subdir/en.dat
      B) basedir/point_subdir/te_*/tau_*/en.dat

    返回：
      layout_mode: "A" or "B"
      points: [(te,tau,dirname),...]  dirname 是 te_.._tau_.. 的真实目录名（在对应 root 下）
    """
    # try A: points are under basedir
    ptsA = scan_te_tau_dirs(basedir)
    if ptsA:
        # 检查第一个点是否存在叶子文件
        te0, tau0, dn0 = ptsA[0]
        leaf0 = os.path.join(basedir, dn0, point_subdir)
        okA = (os.path.exists(os.path.join(leaf0, en_file)) and
               os.path.exists(os.path.join(leaf0, flux_file)))
        if okA:
            return "A", ptsA

    # try B: points are under basedir/point_subdir
    rootB = os.path.join(basedir, point_subdir)
    ptsB = scan_te_tau_dirs(rootB)
    if ptsB:
        te0, tau0, dn0 = ptsB[0]
        leaf0 = os.path.join(rootB, dn0)
        okB = (os.path.exists(os.path.join(leaf0, en_file)) and
               os.path.exists(os.path.join(leaf0, flux_file)))
        if okB:
            return "B", ptsB

    # 如果 A/B 都不行，给出更明确的错误信息
    raise RuntimeError(
        "Cannot detect directory layout.\n"
        f"Checked:\n"
        f"  A) {basedir}/te_*/tau_*/{point_subdir}/{en_file}\n"
        f"  B) {basedir}/{point_subdir}/te_*/tau_*/{en_file}\n"
        "Please verify your basedir and point_subdir."
    )

def point_leaf_dir(basedir: str, layout: str, dirname: str, point_subdir: str) -> str:
    """
    返回包含 en.dat/flux.dat 的叶子目录
    layout A: basedir/dirname/point_subdir
    layout B: basedir/point_subdir/dirname
    """
    if layout == "A":
        return os.path.join(basedir, dirname, point_subdir)
    if layout == "B":
        return os.path.join(basedir, point_subdir, dirname)
    raise ValueError("Unknown layout: " + str(layout))

def load_spectrum_from_leaf(leaf_dir: str, en_file: str, flux_file: str) -> Tuple[np.ndarray, np.ndarray]:
    """从叶子目录读取 (E,F)，横坐标严格使用 en.dat"""
    en_path   = os.path.join(leaf_dir, en_file)
    flux_path = os.path.join(leaf_dir, flux_file)

    if not os.path.exists(en_path):
        raise FileNotFoundError(f"Energy file not found: {en_path}")
    if not os.path.exists(flux_path):
        raise FileNotFoundError(f"Flux file not found: {flux_path}")

    E = load_dat_file(en_path)
    F = load_dat_file(flux_path)

    if E.shape != F.shape:
        raise ValueError(f"Energy/Flux length mismatch in {leaf_dir}: {E.shape} vs {F.shape}")

    if not np.all(np.diff(E) > 0):
        raise ValueError(f"Energy grid not strictly increasing: {leaf_dir}")

    return E, F


# =====================
# math helpers
# =====================
def safe_log10(y):
    y = np.asarray(y, float)
    ok = np.isfinite(y) & (y > 0)
    out = np.full_like(y, np.nan, dtype=float)
    out[ok] = np.log10(np.maximum(y[ok], FLOOR_POS))
    return out

def safe_pow10(x):
    x = np.asarray(x, float)
    out = np.full_like(x, np.nan, dtype=float)
    ok = np.isfinite(x)
    out[ok] = 10.0 ** x[ok]
    return out

def rel_err(pred, true):
    pred = np.asarray(pred, float)
    true = np.asarray(true, float)
    ok = np.isfinite(pred) & np.isfinite(true) & (pred > 0) & (true > 0)
    out = np.full_like(true, np.nan, dtype=float)
    out[ok] = np.abs(1.0 - pred[ok] / true[ok])
    return out

def percentile_stats(err_mat):
    v = err_mat[np.isfinite(err_mat)]
    if v.size == 0:
        return {"median": np.nan, "p90": np.nan, "p95": np.nan, "max": np.nan}
    return {
        "median": float(np.percentile(v, 50)),
        "p90":    float(np.percentile(v, 90)),
        "p95":    float(np.percentile(v, 95)),
        "max":    float(np.nanmax(v)),
    }

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


# =====================
# 1D interpolation kernels
# =====================
def interp_1d_two_point(x1, y1, x2, y2, x):
    t = (x - x1) / (x2 - x1)
    return (1.0 - t) * y1 + t * y2

def predict_one_pt(p_grid, F_grid, pt, mode):
    """
    mode:
      - linp_linF
      - linp_logF
      - logp_linF
      - logp_logF
    """
    j = np.searchsorted(p_grid, pt)
    if j <= 0 or j >= len(p_grid):
        return np.full(F_grid.shape[1], np.nan, dtype=float)

    p1, p2 = p_grid[j-1], p_grid[j]
    y1, y2 = F_grid[j-1], F_grid[j]

    if mode.startswith("logp_"):
        if p1 <= 0 or p2 <= 0 or pt <= 0:
            return np.full(F_grid.shape[1], np.nan, dtype=float)
        x1, x2, x = np.log10(p1), np.log10(p2), np.log10(pt)
    else:
        x1, x2, x = p1, p2, pt

    if mode.endswith("_linF"):
        return interp_1d_two_point(x1, y1, x2, y2, x)

    if mode.endswith("_logF"):
        ly1 = safe_log10(y1)
        ly2 = safe_log10(y2)
        lp  = interp_1d_two_point(x1, ly1, x2, ly2, x)
        return safe_pow10(lp)

    raise ValueError("unknown mode: " + str(mode))


# =====================
# plotting
# =====================
def plot_page_all_spectra(pdf, E, p_vals, F, title_line: str, p_label: str):
    fig, ax = plt.subplots(figsize=(9.5, 7.2), constrained_layout=True)
    for i, (p, y) in enumerate(zip(p_vals, F)):
        if i % 2 == 0:
            ax.loglog(E, y, lw=LINE_W_SOLID, alpha=LINE_ALPHA_SOLID, linestyle="-")
        else:
            ax.loglog(E, y, lw=LINE_W_DASH,  alpha=LINE_ALPHA_DASH,  linestyle="--")

    ax.set_xlabel("Energy [keV] (from en.dat)")
    ax.set_ylabel("Flux(E) (from flux.dat)")
    ax.grid(True, which="both", alpha=0.25)
    ax.set_title(
        "All points in slice: solid=grid(even index), dashed=test(odd index)\n" + title_line
    )
    ax.text(
        0.98, 0.02,
        f"{p_label} points: N={len(p_vals)}  (grid={np.sum(np.arange(len(p_vals))%2==0)}, test={np.sum(np.arange(len(p_vals))%2==1)})",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=10
    )
    pdf.savefig(fig)
    plt.close(fig)

def plot_page_error_heatmap(pdf, E_centers, err_mat, mode, stats):
    E_edges = make_edges_from_centers_log(E_centers)
    y_edges = np.arange(err_mat.shape[0] + 1)

    fig, ax = plt.subplots(figsize=(11.0, 6.6), constrained_layout=True)
    show = np.clip(err_mat, 0, ERR_CLIP)
    show = np.ma.masked_invalid(show)

    im = ax.pcolormesh(E_edges, y_edges, show, shading="auto")
    ax.set_xscale("log")
    ax.set_xlabel("Energy [keV] (from en.dat)")
    ax.set_ylabel("test index (odd)")
    ax.set_title(f"Interpolation error: {mode}   (|1 - pred/true|, clipped to {ERR_CLIP:g})")

    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("relative error")

    txt = (
        f"median={stats['median']:.3g}\n"
        f"p90={stats['p90']:.3g}\n"
        f"p95={stats['p95']:.3g}\n"
        f"max={stats['max']:.3g}"
    )
    ax.text(0.98, 0.98, txt, transform=ax.transAxes, ha="right", va="top", fontsize=10)

    pdf.savefig(fig)
    plt.close(fig)


# =====================
# main
# =====================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--basedir", type=str, default=BASEDIR_DEFAULT, help="warmcom _exports base dir")
    parser.add_argument("--outbase", type=str, default=OUTBASE_DEFAULT, help="output base dir")
    parser.add_argument("--outsub",  type=str, default=OUTSUB_DEFAULT,  help="output sub dir name")

    parser.add_argument("--point-subdir", type=str, default=POINT_SUBDIR_DEFAULT,
                        help="extra layer under each point dir, e.g., calspec (default: calspec)")

    # REQUIRED
    parser.add_argument("--mode", type=str, required=True, choices=["te", "tau"],
                        help="slice mode: te means fix tau and vary te; tau means fix te and vary tau")
    parser.add_argument("--tau", type=float, default=None, help="if mode=te: fixed tau")
    parser.add_argument("--te",  type=float, default=None, help="if mode=tau: fixed te")

    parser.add_argument("--atol", type=float, default=ATOL_DEFAULT, help="slice match atol")
    parser.add_argument("--en-file", type=str, default=EN_FILE_DEFAULT, help="energy binary dat filename")
    parser.add_argument("--flux-file", type=str, default=FLUX_FILE_DEFAULT, help="flux binary dat filename")

    args = parser.parse_args()

    if args.mode == "te" and args.tau is None:
        raise SystemExit("ERROR: --mode te requires --tau <value>")
    if args.mode == "tau" and args.te is None:
        raise SystemExit("ERROR: --mode tau requires --te <value>")

    fixed_tau = float(args.tau) if args.mode == "te" else None
    fixed_te  = float(args.te)  if args.mode == "tau" else None
    atol = float(args.atol)

    # auto-detect layout + scan points (use real dirnames)
    layout, points = detect_layout(args.basedir, args.point_subdir, args.en_file, args.flux_file)

    # slice
    if args.mode == "te":
        slice_pts = [(te, tau, dn) for (te, tau, dn) in points if np.isclose(tau, fixed_tau, rtol=0.0, atol=atol)]
        if len(slice_pts) < 5:
            taus = np.array(sorted({tau for _, tau, _ in points}), float)
            raise RuntimeError(
                f"slice(tau={fixed_tau:g}) has too few points: {len(slice_pts)} (<5)\n"
                f"Available tau (unique) examples: {taus[:20]}{' ...' if len(taus)>20 else ''}\n"
                f"(layout detected: {layout})"
            )
        slice_pts.sort(key=lambda x: x[0])  # by te
        p_vals = np.array([te for te, _, _ in slice_pts], float)
        p_label = "te"
        title_line = f"layout={layout} ; fix tau={fixed_tau:g} vary te ; basedir={args.basedir} ; subdir={args.point_subdir}"
        tag = f"tau{fixed_tau:g}".replace(".", "p")
    else:
        slice_pts = [(te, tau, dn) for (te, tau, dn) in points if np.isclose(te, fixed_te, rtol=0.0, atol=atol)]
        if len(slice_pts) < 5:
            tes = np.array(sorted({te for te, _, _ in points}), float)
            raise RuntimeError(
                f"slice(te={fixed_te:g}) has too few points: {len(slice_pts)} (<5)\n"
                f"Available te (unique) examples: {tes[:20]}{' ...' if len(tes)>20 else ''}\n"
                f"(layout detected: {layout})"
            )
        slice_pts.sort(key=lambda x: x[1])  # by tau
        p_vals = np.array([tau for _, tau, _ in slice_pts], float)
        p_label = "tau"
        title_line = f"layout={layout} ; fix te={fixed_te:g} vary tau ; basedir={args.basedir} ; subdir={args.point_subdir}"
        tag = f"te{fixed_te:g}".replace(".", "p")

    # load spectra for all points (order = p_vals order)
    E_ref = None
    F_list = []
    leaf_dirs = []

    for te, tau, dn in slice_pts:
        leaf = point_leaf_dir(args.basedir, layout, dn, args.point_subdir)
        E, F = load_spectrum_from_leaf(leaf, args.en_file, args.flux_file)

        if E_ref is None:
            E_ref = E.copy()
        else:
            if len(E) != len(E_ref) or (not np.allclose(E, E_ref, rtol=0.0, atol=0.0)):
                raise RuntimeError(f"Energy grid mismatch between points: {leaf}")

        F_list.append(F)
        leaf_dirs.append(leaf)

    F_mat = np.vstack(F_list)

    # split even/odd indices
    idx = np.arange(len(p_vals))
    idx_grid = (idx % 2 == 0)
    idx_test = ~idx_grid

    p_grid, F_grid = p_vals[idx_grid], F_mat[idx_grid]
    p_test, F_test = p_vals[idx_test], F_mat[idx_test]

    # benchmark
    modes = ["linp_linF", "linp_logF", "logp_linF", "logp_logF"]
    err_mats = {}
    stats_all = {}

    for mode in modes:
        errs = []
        for pt, Ft in zip(p_test, F_test):
            P = predict_one_pt(p_grid, F_grid, pt, mode=mode)
            errs.append(rel_err(P, Ft))
        err = np.vstack(errs)
        err_mats[mode] = err
        stats_all[mode] = percentile_stats(err)

    # outputs
    outdir = os.path.join(args.outbase, args.outsub)
    ensure_dir(outdir)

    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_pdf  = os.path.join(outdir, f"warmcom_interp_{args.mode}_{tag}_{ts}.pdf")
    out_txt  = os.path.join(outdir, f"warmcom_interp_{args.mode}_{tag}_{ts}.summary.txt")
    out_used = os.path.join(outdir, f"warmcom_interp_{args.mode}_{tag}_{ts}.used_dirs.txt")

    with open(out_used, "w", encoding="utf-8") as f:
        f.write("# order follows p_vals sorting; leaf_dir  te  tau  dirname  layout\n")
        for (te, tau, dn), leaf in zip(slice_pts, leaf_dirs):
            f.write(f"{leaf}  {te:.16g}  {tau:.16g}  {dn}  {layout}\n")

    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"basedir: {args.basedir}\n")
        f.write(f"layout: {layout}\n")
        f.write(f"point_subdir: {args.point_subdir}\n")
        f.write(f"mode: {args.mode}\n")
        if args.mode == "te":
            f.write(f"fixed_tau: {fixed_tau:.16g}\n")
        else:
            f.write(f"fixed_te: {fixed_te:.16g}\n")
        f.write(f"atol: {atol:.3g}\n")
        f.write(f"en_file: {args.en_file}\n")
        f.write(f"flux_file: {args.flux_file}\n")
        f.write(f"N(all)={len(p_vals)}; grid={len(p_grid)}; test={len(p_test)}\n")
        f.write(f"ERR_CLIP={ERR_CLIP:g}\n")
        f.write(f"used_dirs_list: {out_used}\n\n")
        for mode in modes:
            s = stats_all[mode]
            f.write(f"[{mode}]\n")
            f.write(f"median = {s['median']:.6g}\n")
            f.write(f"p90    = {s['p90']:.6g}\n")
            f.write(f"p95    = {s['p95']:.6g}\n")
            f.write(f"max    = {s['max']:.6g}\n\n")

    with PdfPages(out_pdf) as pdf:
        plot_page_all_spectra(pdf, E_ref, p_vals, F_mat, title_line=title_line, p_label=p_label)
        for mode in modes:
            plot_page_error_heatmap(pdf, E_ref, err_mats[mode], mode, stats_all[mode])

    print("✅ Saved PDF ->", out_pdf)
    print("✅ Saved TXT ->", out_txt)
    print("✅ Saved used-dirs ->", out_used)


if __name__ == "__main__":
    main()
