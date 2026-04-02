#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compare_2tests_warmcom_vs_full.py

一次生成 3 页 PDF：tetest / tautest 相对 FULL 的对比（warmcom slab）。

Pages 1 : tetest  vs FULL   (REF=FULL)
Pages 2 : tautest vs FULL   (REF=FULL)
Pages 3 : FULL    vs FULL   (REF=FULL, sanity check)

输出：
- 3页 PDF
- 1个 txt：unused_params_summary.txt（te/tau 的 excluded 列表 + 最终选取点）

选点规则（模仿你 gdiskwien 那个“excluded 交集”逻辑）：
- 从 WARMCOM_BASEDIR/chose.log 解析所有存在的 (te,tau) 点
- te_all, tau_all 唯一排序
- excluded_te  = te_all[1::2]
- excluded_tau = tau_all[1::2]
- 最终使用 (te,tau) 必须同时在 excluded_te 与 excluded_tau，且该 pair 真实存在于 chose.log

支持：
  --seed
  --te   手动指定 te（必须属于 excluded_te）
  --tau  手动指定 tau（必须属于 excluded_tau）

XSPEC 模型：
  包名/目录：
    warmcom_slab
    warmcom_slab_tetest
    warmcom_slab_tautest
  mo 名称：
    warmcomslab
    warmcomslabtetest
    warmcomslabtautest

注意：
- additive local model 会有 XSPEC 自动 norm 参数，因此 xcm 里要给 4 个参数：
    te, tau, z, norm
- 你的 C 代码只读 (te,tau,z)，norm 会被忽略（但 XSPEC 仍会问）
"""

import os
import re
import uuid
import shutil
import argparse
import subprocess
import numpy as np
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# =========================
# Paths (YOU)
# =========================
OUT_DIR = os.path.expanduser("~/data/monk/plot/warmcorona/test/2026.01.26")

MODEL_BASEDIR = os.path.expanduser("~/data/monk/plot/warmcorona/model/warmcom")
MODEL_DIR_FULL      = os.path.join(MODEL_BASEDIR, "warmcom_slab")
MODEL_DIR_TETEST    = os.path.join(MODEL_BASEDIR, "warmcom_slab_tetest")
MODEL_DIR_TAUTEST   = os.path.join(MODEL_BASEDIR, "warmcom_slab_tautest")

PKG_NAME_FULL     = "warmcom_slab"
MODEL_NAME_FULL   = "warmcomslab"

PKG_NAME_TETEST   = "warmcom_slab_tetest"
MODEL_NAME_TETEST = "warmcomslabtetest"

PKG_NAME_TAUTEST  = "warmcom_slab_tautest"
MODEL_NAME_TAUTEST= "warmcomslabtautest"


# =========================
# Warmcom grid base dir (chose.log lives here)
# =========================
WARMCOM_BASEDIR_DEFAULT = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"
CHOSE_LOG_SUBPATH = "chose.log"


# =========================
# XSPEC energy grid
# =========================
E_MIN = 0.001
E_MAX = 500.0
NBINS = 1000
E_GRID_MODE = "log"

# params
Z_DEFAULT = 0.0
NORM_DEFAULT = 1.0
FLUX_MIN = 1e-80


# -------------------------
# helpers
# -------------------------
def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def skip1(arr, start_idx=1):
    arr = np.asarray(arr, float)
    return arr[start_idx::2]

def parse_chose_log_pairs(basedir: str):
    logpath = os.path.join(basedir, CHOSE_LOG_SUBPATH)
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
    te_all = np.array(sorted({p[0] for p in pairs}), dtype=float)
    tau_all = np.array(sorted({p[1] for p in pairs}), dtype=float)
    return pairs, te_all, tau_all

def write_unused_params_summary_warmcom(exc_te, exc_tau, chosen, out_txt: str):
    ensure_dir(os.path.dirname(out_txt))
    te_use, tau_use = chosen
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("# Unused parameter grids for warmcom tests (derived from chose.log)\n")
        f.write("# Rule: sort unique grid then take odd indices: x[1], x[3], ...\n\n")

        f.write("[warmcomslabtetest]\n")
        f.write("# skipped axis: te\n")
        f.write("te:\n")
        for x in np.asarray(exc_te, float):
            f.write(f"{x:.16g}\n")
        f.write("\n")

        f.write("[warmcomslabtautest]\n")
        f.write("# skipped axis: tau\n")
        f.write("tau:\n")
        for x in np.asarray(exc_tau, float):
            f.write(f"{x:.16g}\n")
        f.write("\n")

        f.write("[chosen point]\n")
        f.write(f"te  = {te_use:.16g}\n")
        f.write(f"tau = {tau_use:.16g}\n")

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

def compute_error_abs(y_model: np.ndarray, y_ref: np.ndarray) -> np.ndarray:
    y_model = np.asarray(y_model, float)
    y_ref   = np.asarray(y_ref, float)
    ok = np.isfinite(y_model) & np.isfinite(y_ref) & (y_model > 0) & (y_ref > 0)
    err = np.full_like(y_model, np.nan, dtype=float)
    err[ok] = np.abs(1.0 - (y_model[ok] / y_ref[ok]))
    return err

def xspec_write_qdp_warmcom(
    model_dir: str, pkg: str, model_name: str,
    te: float, tau: float, z: float, norm: float,
    out_qdp: str, out_xcm: str,
    emin: float, emax: float, nbins: int, grid_mode: str,
    warmcom_basedir: str
):
    ensure_dir(os.path.dirname(out_qdp))
    ensure_dir(os.path.dirname(out_xcm))

    out_dir = os.path.dirname(os.path.abspath(out_qdp))
    tmp_short = f"qdp_tmp_{os.getpid()}_{uuid.uuid4().hex[:6]}.qdp"
    tmp_path  = os.path.join(out_dir, tmp_short)

    grid_kw = "log" if grid_mode.lower().startswith("log") else "lin"

    # IMPORTANT: use iplot -> wdata -> quit, avoid "setplot command exit" killing XSPEC
    cmds = f"""query yes
cpd /null
lmod {pkg} {model_dir}
model {model_name}
{te}
{tau}
{z}
{norm}
setplot rebin 1 1
setplot energy
setplot area off
energies {emin} {emax} {nbins} {grid_kw}
plot model
iplot
wdata {tmp_short}
quit
exit
"""

    with open(out_xcm, "w", encoding="utf-8") as f:
        f.write(cmds)

    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env.setdefault("PGPLOT_DEV", "/null")
    env["WARMCOM_BASEDIR"] = warmcom_basedir  # critical

    proc = subprocess.run(
        ["xspec"],
        input=cmds,
        text=True,
        cwd=out_dir,
        env=env,
        capture_output=True
    )

    if proc.returncode != 0:
        print("\n[XSPEC STDOUT tail]\n", proc.stdout[-3000:])
        print("\n[XSPEC STDERR tail]\n", proc.stderr[-3000:])
        raise RuntimeError(f"XSPEC failed: returncode={proc.returncode}")

    if (not os.path.exists(tmp_path)) or os.path.getsize(tmp_path) == 0:
        print("\n[XSPEC STDOUT tail]\n", proc.stdout[-3000:])
        print("\n[XSPEC STDERR tail]\n", proc.stderr[-3000:])
        raise RuntimeError(f"XSPEC did not produce QDP. Expected: {tmp_path}")

    shutil.move(tmp_path, out_qdp)
    return out_qdp


# -------------------------
# Plotting (same 3-panels style as you had)
# -------------------------
def plot_page(pdf, E, y_model, y_ref,
              label_model, label_ref,
              title, ref_note):
    fig, axes = plt.subplots(
        3, 1, figsize=(8.8, 9.8), sharex=False,
        gridspec_kw={"height_ratios": [2.2, 1.2, 1.4]}
    )
    ax1, ax2, ax3 = axes

    y_model = np.asarray(y_model, float)
    y_ref   = np.asarray(y_ref, float)
    E       = np.asarray(E, float)

    # (1) Flux
    ax1.loglog(E, np.clip(y_model, FLUX_MIN, np.inf), lw=1.8, label=label_model)
    ax1.loglog(E, np.clip(y_ref,   FLUX_MIN, np.inf), ls="--", lw=1.3, label=label_ref)
    ax1.set_ylabel("Flux")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=9)
    ax1.set_ylim(bottom=FLUX_MIN)

    # (2) Abs error vs E
    err_abs = compute_error_abs(y_model, y_ref)
    if np.any(np.isfinite(err_abs)):
        ax2.semilogx(E, err_abs, lw=1.2)
    ax2.set_xlabel("Energy [keV]")
    ax2.set_ylabel(r"$\left|1-\frac{\mathrm{MODEL}}{\mathrm{REF}}\right|$")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.axhline(0.0, ls="--", lw=0.8, color="gray")

    # (3) scatter: x = MODEL/max(MODEL) , y = abs err
    m_ok = np.isfinite(y_model) & (y_model > 0)
    mmax = float(np.nanmax(y_model[m_ok])) if np.any(m_ok) else np.nan
    ok_sc = np.isfinite(err_abs) & np.isfinite(y_model) & (y_model > 0) & np.isfinite(mmax) & (mmax > 0)

    if np.any(ok_sc):
        xratio = y_model[ok_sc] / mmax
        xratio = np.clip(xratio, 1e-12, 1.0)
        ax3.scatter(xratio, err_abs[ok_sc], s=8, alpha=0.75, linewidths=0.0)
        ax3.set_xscale("log")
        ax3.set_xlim(1e-6, 1.05)
    else:
        ax3.text(0.5, 0.5, "No valid points for ratio-error plot",
                 ha="center", va="center", transform=ax3.transAxes)

    ax3.set_xlabel(r"$F_{\rm MODEL}(E)\ /\ \max_E F_{\rm MODEL}(E)$")
    ax3.set_ylabel(r"$\left|1-\frac{\mathrm{MODEL}}{\mathrm{REF}}\right|$")
    ax3.grid(True, which="both", alpha=0.25)

    fig.suptitle(title + "\n" + ref_note, fontsize=11)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# -------------------------
# choose common excluded (te,tau)
# -------------------------
def _in_list(x, arr, atol=1e-12):
    return bool(np.any(np.isclose(np.asarray(arr, float), float(x), rtol=0.0, atol=atol)))

def choose_common_excluded_pair(pairs, exc_te, exc_tau, te_req=None, tau_req=None, seed=None):
    if te_req is not None and (not _in_list(te_req, exc_te)):
        raise RuntimeError(f"你指定的 te={te_req:g} 不在 excluded_te（tetest 未使用列表）里。")
    if tau_req is not None and (not _in_list(tau_req, exc_tau)):
        raise RuntimeError(f"你指定的 tau={tau_req:g} 不在 excluded_tau（tautest 未使用列表）里。")

    cand = []
    for te, tau in pairs:
        if (te_req is not None) and (not np.isclose(te, te_req, rtol=0.0, atol=1e-12)):
            continue
        if (tau_req is not None) and (not np.isclose(tau, tau_req, rtol=0.0, atol=1e-12)):
            continue
        if _in_list(te, exc_te) and _in_list(tau, exc_tau):
            cand.append((te, tau))

    if not cand:
        raise RuntimeError(
            "找不到满足 excluded_te ∩ excluded_tau 的真实 (te,tau) 点。\n"
            "说明：你的 grid 可能不是 te×tau 的完整笛卡尔积，或者 odd-index 的交集刚好为空。\n"
            "解决：你可以只手动指定一维（--te 或 --tau）看看能不能找到交集点。"
        )

    rng = np.random.default_rng(seed)
    return cand[int(rng.integers(0, len(cand)))]


# -------------------------
# main
# -------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=None, help="随机种子（在候选点里随机挑一对）")
    parser.add_argument("--basedir", type=str, default=WARMCOM_BASEDIR_DEFAULT, help="WARMCOM_BASEDIR (contains chose.log)")
    parser.add_argument("--te", type=float, default=None, help="手动指定 te（必须属于 excluded_te）")
    parser.add_argument("--tau", type=float, default=None, help="手动指定 tau（必须属于 excluded_tau）")
    parser.add_argument("--z", type=float, default=Z_DEFAULT, help="redshift z")
    parser.add_argument("--norm", type=float, default=NORM_DEFAULT, help="XSPEC norm (model ignores it)")
    args = parser.parse_args()

    ensure_dir(OUT_DIR)

    pairs, te_all, tau_all = parse_chose_log_pairs(args.basedir)
    exc_te  = skip1(te_all,  1)
    exc_tau = skip1(tau_all, 1)

    chosen_te, chosen_tau = choose_common_excluded_pair(
        pairs, exc_te, exc_tau,
        te_req=args.te, tau_req=args.tau,
        seed=args.seed
    )

    unused_txt = os.path.join(OUT_DIR, "unused_params_summary.txt")
    write_unused_params_summary_warmcom(exc_te, exc_tau, (chosen_te, chosen_tau), unused_txt)
    print("Unused-parameter summary written ->", unused_txt)

    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    tag = f"warmcom_compare_2tests_onePoint_te{chosen_te:g}_tau{chosen_tau:g}_s{args.seed if args.seed is not None else 'none'}_{ts}"
    out_pdf = os.path.join(OUT_DIR, f"{tag}.pdf")

    print("\n== ONE COMMON PARAM SET (used for ALL pages) ==")
    print(f"  te   = {chosen_te:g}   (excluded_te)")
    print(f"  tau  = {chosen_tau:g}  (excluded_tau)")
    print(f"  z    = {args.z:g}, norm={args.norm:g}")
    print(f"  basedir = {args.basedir}")
    print(f"  PDF OUT = {out_pdf}")

    # write QDPs
    qdp_full   = os.path.join(OUT_DIR, f"{tag}__p03_full.qdp")
    xcm_full   = os.path.join(OUT_DIR, f"{tag}__p03_full.xcm")

    qdp_tetest = os.path.join(OUT_DIR, f"{tag}__p01_tetest.qdp")
    xcm_tetest = os.path.join(OUT_DIR, f"{tag}__p01_tetest.xcm")

    qdp_tautest= os.path.join(OUT_DIR, f"{tag}__p02_tautest.qdp")
    xcm_tautest= os.path.join(OUT_DIR, f"{tag}__p02_tautest.xcm")

    # FULL
    xspec_write_qdp_warmcom(
        MODEL_DIR_FULL, PKG_NAME_FULL, MODEL_NAME_FULL,
        chosen_te, chosen_tau, float(args.z), float(args.norm),
        qdp_full, xcm_full,
        E_MIN, E_MAX, NBINS, E_GRID_MODE,
        warmcom_basedir=args.basedir
    )
    Efull, Ffull = read_qdp_model(qdp_full)

    # TETEST
    xspec_write_qdp_warmcom(
        MODEL_DIR_TETEST, PKG_NAME_TETEST, MODEL_NAME_TETEST,
        chosen_te, chosen_tau, float(args.z), float(args.norm),
        qdp_tetest, xcm_tetest,
        E_MIN, E_MAX, NBINS, E_GRID_MODE,
        warmcom_basedir=args.basedir
    )
    Et, Ft = read_qdp_model(qdp_tetest)

    # TAUTEST
    xspec_write_qdp_warmcom(
        MODEL_DIR_TAUTEST, PKG_NAME_TAUTEST, MODEL_NAME_TAUTEST,
        chosen_te, chosen_tau, float(args.z), float(args.norm),
        qdp_tautest, xcm_tautest,
        E_MIN, E_MAX, NBINS, E_GRID_MODE,
        warmcom_basedir=args.basedir
    )
    Eu, Fu = read_qdp_model(qdp_tautest)

    # align onto FULL energy grid
    Ft_on = np.interp(Efull, Et, Ft, left=np.nan, right=np.nan)
    Fu_on = np.interp(Efull, Eu, Fu, left=np.nan, right=np.nan)

    with PdfPages(out_pdf) as pdf:
        # Page 1: tetest vs FULL
        plot_page(
            pdf, Efull, Ft_on, Ffull,
            "warmcomslabtetest (XSPEC)", "warmcomslab FULL (XSPEC)",
            title=(f"Page 1 / tetest vs FULL\n"
                   f"te={chosen_te:g}; tau={chosen_tau:g}; z={args.z:g}; norm={args.norm:g}"),
            ref_note="REF = FULL"
        )

        # Page 2: tautest vs FULL
        plot_page(
            pdf, Efull, Fu_on, Ffull,
            "warmcomslabtautest (XSPEC)", "warmcomslab FULL (XSPEC)",
            title=(f"Page 2 / tautest vs FULL\n"
                   f"te={chosen_te:g}; tau={chosen_tau:g}; z={args.z:g}; norm={args.norm:g}"),
            ref_note="REF = FULL"
        )

        # Page 3: FULL vs FULL (sanity)
        plot_page(
            pdf, Efull, Ffull, Ffull,
            "warmcomslab FULL (XSPEC)", "warmcomslab FULL (XSPEC)",
            title=(f"Page 3 / FULL vs FULL (sanity check)\n"
                   f"te={chosen_te:g}; tau={chosen_tau:g}; z={args.z:g}; norm={args.norm:g}"),
            ref_note="REF = FULL"
        )

    print("\n✅ DONE")
    print("  PDF =", out_pdf)
    print("  unused summary txt =", unused_txt)


if __name__ == "__main__":
    main()
