#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare nthComp vs THREE warmcom models (same Te, tau) and overlay
best-fit Gamma straight lines (log–log fit) on the FIRST page.

Enhancements:
- Ensure slab is included: assert all .mod files exist before running.
- Dynamic lower bound for fit: per spectrum starts from first valid energy point.
- Upper bound fixed to 1.0 keV.
- Slab follows the exact same pipeline as spheres (normalize, fit, overlay, residuals).
- EXTRA: At the end of the same PDF, replot already-computed XSPEC
  model curves:
    * Page: nthComp vs w1_sphere_107
    * Page: nthComp vs w3_slab_106
  using exactly the same data arrays and axes as page 1
  (no new XSPEC runs, no新模型).

Keeps the user's original nthComp XCM logic (with newpar):
model nthComp
{gamma}
{max(kTe_nth, 0.01)}
{kTbb_nth}
0
{z_value}
1
newpar 2 {kTe_nth} 1e-3 0.1 0.1 10 1000
"""

import os, sys, subprocess, numpy as np
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from shutil import which

# ==================== User config ====================
TARGET_DIR = "/home/hdw/data/monk/plot/warmcorona/test/compare_warmcom_nth"

# Three warmcom models (paths you provided) — slab INCLUDED
WARMCOM_MODELS = [
    ("w1_sphere_107",
     "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-1.0_5-25_smoothed_tv_Tom_new_107.mod",
     "sphere"),
    ("w2_sphere_106",
     "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-2.0_5-22_smoothed_tv_clean_Tom_new.mod",
     "sphere"),
    ("w3_slab_106",
     "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_pure_Tom_slab_107.mod",
     "slab"),
]

# Normalization and fit settings
E_REF = 0.15                      # normalization energy (keV)
E_MAX_FIT = 1.0                   # Gamma fit upper bound; lower bound is dynamic per spectrum
E_MIN, E_MAX, NPTS = 0.1, 10.0, 1000

# ==================== Inputs ====================
def input_with_default(prompt, default_str):
    try:
        s = input(f"{prompt} [default {default_str}]: ").strip()
        return s if s else default_str
    except EOFError:
        return default_str

def as_float(s):
    try:
        return float(s)
    except Exception:
        raise RuntimeError(f"Need a numeric value: {s}")

kTe_nth   = as_float(input_with_default("nthComp kT_e (keV)", "0.5"))
gamma     = as_float(input_with_default("nthComp Gamma", "2.5"))
z_value   = as_float(input_with_default("redshift z", "0.0"))
kTbb_nth  = as_float(input_with_default("nthComp kTbb (keV)", "0.003"))

# ==================== Preflight: ensure all .mod exist (including slab) ====================
missing = [p for _, p, _ in WARMCOM_MODELS if not os.path.exists(p)]
if missing:
    raise FileNotFoundError(
        "These .mod files are missing (slab must exist too):\n" +
        "\n".join(missing)
    )

# ==================== Derive tau from (Gamma, kTe) ====================
# 使用 Pozdnyakov et al. (1977) 形式
# theta_e = kTe / 511 keV
# term_g  = (Gamma + 0.5)^2 - 9/4
# tau     = sqrt( 9/4 + 3 / (theta_e * term_g) ) - 3/2
me_keV  = 511.0
theta_e = kTe_nth / me_keV
term_g  = (gamma + 0.5)**2 - 2.25  # 9/4

if theta_e <= 0 or term_g <= 0:
    raise RuntimeError("Unphysical parameters: need kT_e>0 and (Gamma+0.5)^2 > 9/4.")

tau_nth = (2.25 + 3.0 / (theta_e * term_g))**0.5 - 1.5
te_warm, tau_warm = kTe_nth, tau_nth

print(f"[INFO] Derived tau from (Gamma={gamma:.4f}, kTe={kTe_nth:.4f} keV): tau ≈ {tau_nth:.6f}")

# ==================== Helpers ====================
os.makedirs(TARGET_DIR, exist_ok=True)
have_xspec = which("xspec") is not None

def rm_if_exists(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception as e:
        print(f"[WARN] Cannot remove existing file {path}: {e}")

def write_xcm_nthcomp(tag, gamma, kTe, kTbb, z, out_dir):
    """
    Exactly follow user's original nthComp snippet (keep newpar line),
    plus energies/plot/export.
    """
    xcm = os.path.join(out_dir, f"{tag}.xcm")
    qdp = os.path.join(out_dir, f"{tag}.qdp")
    rm_if_exists(qdp)
    xcm_txt = f"""
model nthComp
{gamma}
{max(kTe, 0.01)}
{kTbb}
0
{z}
1
newpar 2 {kTe} 1e-3 0.1 0.1 10 1000
cpd /null
energies 0.1 10 1000 log
plot model
setplot command wdata {qdp}
setplot command exit
iplot
exit
""".lstrip()
    with open(xcm, "w") as f:
        f.write(xcm_txt)
    return xcm, qdp

def write_xcm_warmcom(tag, mod_path, te, tau, z, out_dir):
    """
    Independent XCM/QDP per warmcom (short filenames; delete QDP first).
    """
    xcm = os.path.join(out_dir, f"{tag}.xcm")
    qdp = os.path.join(out_dir, f"{tag}.qdp")
    rm_if_exists(qdp)
    xcm_txt = f"""
model atable{{{mod_path}}}
{te}
{tau}
{z}
1
energies 0.1 10 1000 log
cpd /null
plot model
setplot command wdata {qdp}
setplot command exit
iplot
exit
""".lstrip()
    with open(xcm, "w") as f:
        f.write(xcm_txt)
    return xcm, qdp

def run_xspec(xcm):
    try:
        subprocess.run(
            ["xspec", "-", xcm],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        return True
    except Exception as e:
        print(f"[WARN] XSPEC failed: {e}")
        return False

def read_qdp(qdp_file):
    x, y = [], []
    if not os.path.exists(qdp_file):
        return np.array([]), np.array([])
    with open(qdp_file, "r") as f:
        for line in f:
            if line.startswith(("!", "@", "READ")):
                continue
            cols = line.split()
            if len(cols) >= 3:
                try:
                    x.append(float(cols[0]))
                    y.append(float(cols[2]))
                except Exception:
                    pass
    return np.array(x, dtype=float), np.array(y, dtype=float)

def interp_safe(x_new, x_old, y_old):
    return np.interp(x_new, x_old, y_old, left=np.nan, right=np.nan)

def first_valid_energy(E, F):
    """Return the first energy where spectrum is valid (>0 and finite)."""
    m = np.isfinite(E) & np.isfinite(F) & (E > 0) & (F > 0)
    if not np.any(m):
        return np.nan
    return float(E[m][0])

def range_fit_loglog(E, F, Emin, Emax, min_pts=10):
    """
    Fit log F vs log E in [Emin, Emax], return (Gamma, m, b),
    where m, b are slope/intercept in log-space: logF = m*logE + b, Gamma = -m.
    """
    mask = (
        np.isfinite(E) & np.isfinite(F) &
        (E > 0) & (F > 0) &
        (E >= Emin) & (E <= Emax)
    )
    E2, F2 = E[mask], F[mask]
    if E2.size < min_pts or Emin >= Emax:
        return np.nan, np.nan, np.nan
    x, y = np.log(E2), np.log(F2)
    m, b = np.polyfit(x, y, 1)
    return -m, m, b

# ==================== Generate XCM/QDP ====================
nth_tag = "n"
xcm_n, qdp_n = write_xcm_nthcomp(nth_tag, gamma, kTe_nth, kTbb_nth, z_value, TARGET_DIR)

warm_jobs = []
for idx, (name, path, kind) in enumerate(WARMCOM_MODELS, start=1):
    tag = f"w{idx}"  # w1/w2/w3
    xcm_w, qdp_w = write_xcm_warmcom(tag, path, te_warm, tau_warm, z_value, TARGET_DIR)
    warm_jobs.append((name, kind, path, tag, xcm_w, qdp_w))

if have_xspec:
    ok = run_xspec(xcm_n)
    if not ok:
        sys.exit("[ERROR] XSPEC failed on nthComp.")
    for nm, kd, pth, tg, xcm_w, _ in warm_jobs:
        if not run_xspec(xcm_w):
            sys.exit(f"[ERROR] XSPEC failed on {nm} ({kd}).")
else:
    print("[HINT] XSPEC not found in PATH; only XCMs are created.")
    sys.exit(0)

# ==================== Read QDP & compute/plot ====================
E_n, F_n = read_qdp(qdp_n)
if E_n.size == 0:
    sys.exit("[ERROR] nthComp QDP not found/readable.")

# Build PDF
F_n_ref = np.interp(E_REF, E_n, F_n, left=np.nan, right=np.nan)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
pdf_name = (
    f"cmp_warm_vs_nth_Te{str(kTe_nth).replace('.','p')}"
    f"_G{str(gamma).replace('.','p')}"
    f"_z{str(z_value).replace('.','p')}_{timestamp}.pdf"
)
pdf_path = os.path.join(TARGET_DIR, pdf_name)

# Prepare results dicts (with fitted Gamma on plotted arrays)
results = []

# nthComp fit (on its native array) with dynamic Emin
Emin_n = first_valid_energy(E_n, F_n)
Gamma_nth, m_nth, b_nth = range_fit_loglog(E_n, F_n, Emin_n, E_MAX_FIT)
print(f"[nthComp] fit range = [{Emin_n:.3f}, {E_MAX_FIT:.3f}] keV | Gamma = {Gamma_nth:.4f}")

results.append({
    "label": "nthComp",
    "color": "C0",
    "E": E_n,
    "F": F_n,
    "Gamma_range": Gamma_nth,
    "m": m_nth,
    "b": b_nth,
    "Emin_fit": Emin_n,
    "scaled": False
})

# Warmcom: read, scale to nthComp at E_REF, resample to E_n for plotting, then fit on plotted array
colors = ["C1", "C2", "C3"]
for i, (name, kind, path, tag, _, qdp_w) in enumerate(warm_jobs):
    E_w, F_w = read_qdp(qdp_w)
    if E_w.size == 0:
        sys.exit(f"[ERROR] Missing QDP for {name} ({kind}); slab must be included.")
    F_w_ref = np.interp(E_REF, E_w, F_w, left=np.nan, right=np.nan)
    if not (np.isfinite(F_n_ref) and np.isfinite(F_w_ref) and F_w_ref != 0):
        sys.exit(f"[ERROR] Bad normalization at E_ref for {name} ({kind}).")
    scale = F_n_ref / F_w_ref

    # Resample scaled model onto E_n for consistent plotting & fitting
    F_w_scaled_on_En = interp_safe(E_n, E_w, scale * F_w)

    Emin_w = first_valid_energy(E_n, F_w_scaled_on_En)
    Gamma_w, m_w, b_w = range_fit_loglog(E_n, F_w_scaled_on_En, Emin_w, E_MAX_FIT)
    print(
        f"[{name}] kind={kind} | scale={scale:.4g} | "
        f"fit range=[{Emin_w:.3f},{E_MAX_FIT:.3f}] keV | Gamma = {Gamma_w:.4f}"
    )

    results.append({
        "label": f"{name} ({kind}) × {scale:.3g}",
        "color": colors[i % len(colors)],
        "E": E_n,  # already on E_n grid
        "F": F_w_scaled_on_En,
        "Gamma_range": Gamma_w,
        "m": m_w,
        "b": b_w,
        "Emin_fit": Emin_w,
        "scaled": True
    })

# ====== PDF generation ======
XLIM_MAIN = None  # x-limits from first page
YLIM_MAIN = None  # y-limits from first page

param_text = (
    f"nthComp: kTe={kTe_nth:.3f} keV, Gamma(input)={gamma:.3f}, "
    f"kTbb={kTbb_nth:.4f}, z={z_value}\n"
    f"tau(derived)≈{tau_nth:.4f}, normalize@{E_REF:.3f} keV\n"
    f"Dashed lines: per-spectrum log–log fits over [Emin_data, {E_MAX_FIT:.1f}] keV."
)

with PdfPages(pdf_path) as pdf:
    # ====== Page 1: Spectra + fitted straight lines ======
    fig1 = plt.figure()
    # Plot spectra
    for r in results:
        lw = 1.8 if r["label"] == "nthComp" else 1.2
        plt.loglog(
            r["E"], r["F"],
            color=r["color"], lw=lw,
            label=f"{r['label']} | Gamma_fit={r['Gamma_range']:.3f}"
        )
    # Overlay fitted straight lines over [Emin_fit, E_MAX_FIT]
    for r in results:
        if (
            not np.isfinite(r["m"]) or
            not np.isfinite(r["b"]) or
            not np.isfinite(r["Emin_fit"])
        ):
            continue
        if r["Emin_fit"] >= E_MAX_FIT:
            continue
        Efit = np.geomspace(r["Emin_fit"], E_MAX_FIT, 200)
        yfit = np.exp(r["m"] * np.log(Efit) + r["b"])  # linear for plotting
        plt.loglog(Efit, yfit, ls="--", lw=1.1, color=r["color"])

    plt.axvline(E_REF, ls="--", lw=0.8, color="gray")
    plt.xlabel("Energy [keV]")
    plt.ylabel("Model spectrum")
    plt.legend(fontsize=7)
    plt.title(
        f"warmcom vs nthComp (same Te, tau); dashed = fit over [Emin_data, {E_MAX_FIT:.1f}] keV"
    )

    # Record axes limits from the first page for later use
    XLIM_MAIN = plt.xlim()
    YLIM_MAIN = plt.ylim()

    fig1.text(
        0.01, 0.98, param_text, ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85),
        fontsize=8
    )
    pdf.savefig(fig1)
    plt.close(fig1)

    # ====== Page 2+: Residuals per model ======
    nth_E, nth_F = results[0]["E"], results[0]["F"]
    for r in results[1:]:
        residual = nth_F - r["F"]
        ratio    = r["F"] / nth_F
        fig, ax1 = plt.subplots()
        ax1.semilogx(nth_E, residual, label="Residual (nth - model×s)")
        ax1.axhline(0, ls="--", lw=0.8)
        ax2 = ax1.twinx()
        ax2.semilogx(nth_E, ratio, 'r', alpha=0.75, label="Ratio (model/nth)")
        ax2.axhline(1.0, ls=":", lw=0.8)
        ax1.set_xlabel("Energy [keV]")
        ax1.set_ylabel("Residual")
        ax2.set_ylabel("Ratio")
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, loc="best", fontsize=8)
        info = (
            f"{r['label']} | Gamma_fit[{r['Emin_fit']:.3f}-{E_MAX_FIT:.3f}]="
            f"{r['Gamma_range']:.3f}\n"
            f"{param_text}"
        )
        fig.text(
            0.01, 0.98, info, ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85),
            fontsize=8
        )
        ax1.axvline(E_REF, ls="--", lw=0.8)
        ax2.axvline(E_REF, ls="--", lw=0.8)
        pdf.savefig(fig)
        plt.close(fig)

    # ====== EXTRA: two separate figures reusing already computed data ======
    # 目标：只用 results 中的数组，重新组合画 nth + w1 / nth + w3，
    # 坐标范围和第一页完全一致，不再跑 XSPEC，也不再新建模型。
    target_names = ["w1_sphere_107", "w3_slab_106"]
    for tname in target_names:
        # 找到对应的 warmcom 结果（已经在第一页画过）
        r_model = None
        for r in results[1:]:
            if r["label"].startswith(tname):
                r_model = r
                break
        if r_model is None:
            continue

        fig = plt.figure()

        # 直接重用 results 里的 E/F 数据（前面 XSPEC 已经跑过）
        # nthComp
        plt.loglog(
            results[0]["E"], results[0]["F"],
            color=results[0]["color"], lw=1.8,
            label=f"{results[0]['label']} | Gamma_fit={results[0]['Gamma_range']:.3f}"
        )
        # 对应 warmcom
        plt.loglog(
            r_model["E"], r_model["F"],
            color=r_model["color"], lw=1.4,
            label=f"{r_model['label']} | Gamma_fit={r_model['Gamma_range']:.3f}"
        )

        # 拟合直线也只用前面算好的 m/b，在同一坐标系上画一遍
        for r in (results[0], r_model):
            if (
                np.isfinite(r["m"]) and np.isfinite(r["b"])
                and np.isfinite(r["Emin_fit"]) and r["Emin_fit"] < E_MAX_FIT
            ):
                Efit = np.geomspace(max(0.1, r["Emin_fit"]), E_MAX_FIT, 200)
                yfit = np.exp(r["m"] * np.log(Efit) + r["b"])
                plt.loglog(Efit, yfit, ls="--", lw=1.0, color=r["color"])

        plt.axvline(E_REF, ls="--", lw=0.8, color="gray")
        plt.xlabel("Energy [keV]")
        plt.ylabel("Model spectrum")

        # 坐标范围：严格和第一页一样
        if XLIM_MAIN is not None:
            plt.xlim(XLIM_MAIN)
        if YLIM_MAIN is not None:
            plt.ylim(YLIM_MAIN)

        title_txt = f"nthComp vs {tname} (same XSPEC curves & axes as page 1)"
        plt.title(title_txt)

        # 下边留出空间写文字
        fig.subplots_adjust(bottom=0.28)

        info_txt = (
            f"{tname}\n"
            f"nthComp fit: Gamma={results[0]['Gamma_range']:.3f} "
            f"[{results[0]['Emin_fit']:.3f}-{E_MAX_FIT:.3f}] keV\n"
            f"{tname} fit: Gamma={r_model['Gamma_range']:.3f} "
            f"[{r_model['Emin_fit']:.3f}-{E_MAX_FIT:.3f}] keV\n"
            f"{param_text}"
        )
        # 文本在图下方，不遮挡曲线
        fig.text(
            0.01, 0.03, info_txt, ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9),
            fontsize=8
        )

        plt.legend(fontsize=7)
        pdf.savefig(fig)
        plt.close(fig)

print(f"[OK] PDF saved: {pdf_path}")
print("[DONE]")
