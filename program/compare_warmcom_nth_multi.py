#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare multiple warmcom models (sphere/slab) to nthComp on one PDF.
- Generates unique XCM + QDP for each model (avoids XSPEC reusing same file).
- Normalizes warmcom to nthComp at E_ref (default 0.15 keV).
- Saves combined PDF with informative filename in the requested folder.

Requirements:
- XSPEC in PATH (optional but needed to actually generate QDPs)
- Provided .mod files exist
"""

import os
import sys
import subprocess
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from shutil import which

# ---------------------- User configuration ----------------------
# Target output directory for XCM/QDP/PDF:
TARGET_DIR = "/home/hdw/data/monk/plot/warmcorona/test/compare_warmcom_nth"

# Warmcom model list: (short_name, path_to_mod, kind)
# You gave two SPHERE models (one "107", one "106"). If you have a SLAB 106 model,
# add it here similarly (uncomment and set the correct path).
WARMCOM_MODELS = [
    ("sphere_107",
     "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-1.0_5-25_smoothed_tv_clean_Tom_new_107.mod",
     "sphere"),
    ("sphere_106",
     "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-2.0_5-22_smoothed_tv_clean_Tom_new.mod",
     "sphere"),
    # Example (edit path to your slab .mod if available):
    # ("slab_106",
    #  "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_slab_0.1-2.0_5-22_smoothed_tv_clean_Tom_new_106.mod",
    #  "slab"),
]

# Reference energy (keV) for normalization:
E_REF = 0.15

# Analysis energy grid for plotting & XSPEC 'energies' command:
E_MIN, E_MAX, NPTS = 0.1, 10.0, 1000  # log grid

# ---------------------- Inputs (interactive) ----------------------
def input_with_default(prompt, default_str):
    try:
        s = input(f"{prompt} [默认 {default_str}]: ").strip()
        return s if s else default_str
    except EOFError:
        return default_str

def as_float(s):
    try:
        return float(s)
    except:
        raise RuntimeError(f"参数需要是数字：{s}")

kTe_nth   = as_float(input_with_default("请输入 nthComp 的 kT_e (keV)", "0.5"))
gamma     = as_float(input_with_default("请输入 nthComp 的 Γ", "2.5"))
z_value   = as_float(input_with_default("请输入 redshift z", "0.0"))
kTbb_nth  = as_float(input_with_default("请输入 nthComp 的 kTbb (keV)", "0.003"))

# ---------------------- Derived tau from (Gamma,kTe) ----------------------
me_keV  = 511.0
theta_e = kTe_nth / me_keV
term_g  = (gamma + 0.5)**2 - 2.25  # 9/4
if theta_e <= 0 or term_g <= 0:
    raise RuntimeError("参数不物理：需满足 kT_e>0 且 (Gamma+0.5)^2 > 9/4。")
tau_nth = (1.0 / (theta_e * term_g))**0.5 - 1.0/3.0

te_warm  = kTe_nth
tau_warm = tau_nth

# ---------------------- Paths & safety ----------------------
os.makedirs(TARGET_DIR, exist_ok=True)
have_xspec = which("xspec") is not None

def log(msg):
    print(msg, flush=True)

def make_log_grid(emin, emax, npts):
    return np.geomspace(emin, emax, npts)

E_grid = make_log_grid(E_MIN, E_MAX, NPTS)

# ---------------------- XSPEC helpers ----------------------
def write_xcm_for_warmcom(short_name, mod_path, te, tau, z, energies, out_dir):
    qdp_file = os.path.join(out_dir, f"{short_name}.qdp")
    xcm_file = os.path.join(out_dir, f"{short_name}.xcm")
    energies_line = f"energies {energies[0]} {energies[-1]} {len(energies)} log"
    xcm = f"""
model atable{{{mod_path}}}
{te}
{tau}
{z}
1
{energies_line}
cpd /null
plot model
setplot command wdata {qdp_file}
setplot command exit
iplot
exit
""".lstrip()
    with open(xcm_file, "w") as f:
        f.write(xcm)
    return xcm_file, qdp_file

def write_xcm_for_nthcomp(short_name, gamma, kTe, kTbb, z, energies, out_dir):
    qdp_file = os.path.join(out_dir, f"{short_name}.qdp")
    xcm_file = os.path.join(out_dir, f"{short_name}.xcm")
    energies_line = f"energies {energies[0]} {energies[-1]} {len(energies)} log"
    xcm = f"""
model nthComp
{gamma}
{max(kTe, 0.01)}
{kTbb}
0
{z}
1
newpar 2 {kTe} 1e-3 0.1 0.1 10 1000
cpd /null
{energies_line}
plot model
setplot command wdata {qdp_file}
setplot command exit
iplot
exit
""".lstrip()
    with open(xcm_file, "w") as f:
        f.write(xcm)
    return xcm_file, qdp_file

def run_xspec_script(xcm_path):
    try:
        subprocess.run(["xspec", "-", xcm_path], check=True)
        return True
    except Exception as e:
        log(f"[警告] 运行 XSPEC 失败：{e}")
        return False

# ---------------------- QDP reading ----------------------
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
                except:
                    pass
    return np.array(x, dtype=float), np.array(y, dtype=float)

# ---------------------- Prepare & run ----------------------
# Write nthComp
nth_tag = f"nthcomp_kTe{str(kTe_nth).replace('.','p')}_G{str(gamma).replace('.','p')}_z{str(z_value).replace('.','p')}"
xcm_nth, qdp_nth = write_xcm_for_nthcomp(nth_tag, gamma, kTe_nth, kTbb_nth, z_value, E_grid, TARGET_DIR)

# Write warmcom XCMs
warm_jobs = []
for short_name, mod_path, kind in WARMCOM_MODELS:
    tag = f"{kind}_{short_name}_te{str(te_warm).replace('.','p')}_tau{str(tau_warm).replace('.','p')}_z{str(z_value).replace('.','p')}"
    xcm_wc, qdp_wc = write_xcm_for_warmcom(tag, mod_path, te_warm, tau_warm, z_value, E_grid, TARGET_DIR)
    warm_jobs.append((short_name, kind, mod_path, xcm_wc, qdp_wc))

# Run XSPEC if available
if have_xspec:
    log("[信息] 检测到 XSPEC，将生成 QDP。")
    run_xspec_script(xcm_nth)
    for _, _, _, xcm_wc, _ in warm_jobs:
        run_xspec_script(xcm_wc)
else:
    log("[提示] 未检测到 XSPEC：将跳过运行，仅生成 XCM 脚本。请在有 XSPEC 的环境运行本脚本以生成 QDP 和 PDF。")

# ---------------------- Load QDP and plot ----------------------
E_nth, F_nth = read_qdp(qdp_nth)
if E_nth.size == 0:
    log("[警告] 没有读取到 nthComp 的 QDP 数据；无法绘图。")
    sys.exit(0)

# Build PDF filename
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
pdf_name = f"compare_warm_vs_nth_kTe{str(kTe_nth).replace('.','p')}_G{str(gamma).replace('.','p')}_z{str(z_value).replace('.','p')}_{timestamp}.pdf"
pdf_path = os.path.join(TARGET_DIR, pdf_name)

def interp_safe(x_new, x_old, y_old):
    return np.interp(x_new, x_old, y_old, left=np.nan, right=np.nan)

# Prepare figure
param_text = (
    f"nthComp: kTe={kTe_nth:.3f} keV, Γ={gamma:.3f}, kTbb={kTbb_nth:.4f} keV, z={z_value}\n"
    f"→ tau(approx) ≈ {tau_nth:.4f} | normalize at {E_REF:.3f} keV\n"
    f"Warmcom models: {', '.join([nm for nm,_,_ in WARMCOM_MODELS])}"
)

with PdfPages(pdf_path) as pdf:
    # Page 1: All model spectra (normalized at E_REF)
    fig1 = plt.figure()
    plt.loglog(E_nth, F_nth, label="nthComp", lw=1.6)
    # normalization reference
    F_n_ref = np.interp(E_REF, E_nth, F_nth, left=np.nan, right=np.nan)

    colors = ["C1","C2","C3","C4","C5","C6"]
    for idx, (short_name, kind, mod_path, xcm_wc, qdp_wc) in enumerate(warm_jobs):
        E_w, F_w = read_qdp(qdp_wc)
        if E_w.size == 0:
            continue
        F_w_ref = np.interp(E_REF, E_w, F_w, left=np.nan, right=np.nan)
        if not (np.isfinite(F_n_ref) and np.isfinite(F_w_ref) and F_w_ref != 0):
            continue
        scale = F_n_ref / F_w_ref
        F_w_scaled_on_nthE = interp_safe(E_nth, E_w, scale*F_w)
        plt.loglog(E_nth, F_w_scaled_on_nthE, color=colors[idx%len(colors)],
                   label=f"{short_name} ({kind}) × {scale:.3g}", lw=1.2)
    plt.axvline(E_REF, ls="--", lw=0.8)
    plt.xlabel("Energy [keV]"); plt.ylabel("Model spectrum")
    plt.legend(loc="best", fontsize=8)
    fig1.text(0.01, 0.98, param_text, ha="left", va="top",
              bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.85), fontsize=8)
    pdf.savefig(fig1); plt.close(fig1)

    # Page 2+: Residuals per model
    for idx, (short_name, kind, mod_path, xcm_wc, qdp_wc) in enumerate(warm_jobs):
        E_w, F_w = read_qdp(qdp_wc)
        if E_w.size == 0:
            continue
        F_w_ref = np.interp(E_REF, E_w, F_w, left=np.nan, right=np.nan)
        if not (np.isfinite(F_n_ref) and np.isfinite(F_w_ref) and F_w_ref != 0):
            continue
        scale = F_n_ref / F_w_ref
        F_w_scaled_on_nthE = interp_safe(E_nth, E_w, scale*F_w)
        if not np.all(np.isfinite(F_w_scaled_on_nthE)):
            continue
        residual = F_nth - F_w_scaled_on_nthE
        ratio    = F_w_scaled_on_nthE / F_nth

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.semilogx(E_nth, residual, label=f"Residual (nth - {short_name}×s)", lw=1.2)
        ax1.axhline(0, ls="--", lw=0.8)
        ax1.set_xlabel("Energy [keV]"); ax1.set_ylabel("Residual")
        ax2 = ax1.twinx()
        ax2.semilogx(E_nth, ratio, label="Ratio (model/nth)", lw=1.0, alpha=0.75)
        ax2.axhline(1.0, ls=":", lw=0.8); ax2.set_ylabel("Ratio")
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines+lines2, labels+labels2, loc="best", fontsize=8)
        fig.text(0.01, 0.98,
                 f"{short_name} ({kind}) | scale @ {E_REF:.3f} keV = {scale:.4g}\n"
                 f"{param_text}",
                 ha="left", va="top",
                 bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.85), fontsize=8)
        ax1.axvline(E_REF, ls="--", lw=0.8); ax2.axvline(E_REF, ls="--", lw=0.8)
        pdf.savefig(fig); plt.close(fig)

print(f"[OK] XCM/QDP/PDF 目标目录：{TARGET_DIR}")
print(f"[OK] 已写出：{os.path.basename(pdf_path)}")
print("[提示] 若未检测到 XSPEC，请在有 XSPEC 的机器上再次运行本脚本以生成 QDP 和 PDF。")
