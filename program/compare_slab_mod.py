#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare THREE warmcom slab models (clean/raw/pure) at the same (Te, tau),
with per-spectrum log–log Gamma fits (dynamic Emin; fixed Emax=1.0 keV).
No nthComp involved.

Pipeline:
- Assert three .mod exist
- For each model: build XCM, XSPEC -> QDP (0.1–10 keV log, 1000 pts)
- Use CLEAN as the reference spectrum
- Normalize RAW and PURE to CLEAN at E_REF
- Fit log F vs log E over [Emin_data, E_MAX_FIT] per spectrum
- PDF:
  * Page 1: 3 spectra + dashed fitted lines
  * Page 2-3: residual (clean - model×s) and ratio (model×s / clean)
"""

import os, sys, subprocess, numpy as np
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from shutil import which

# ==================== User config ====================
TARGET_DIR = "/home/hdw/data/monk/plot/warmcorona/test/compare_warmcom_slab_triplet"

# Three slab models
TABLEFILE_CLEAN = r"/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_smoothed_tv_clean_Tom_slab_107.mod"
TABLEFILE_RAW   = r"/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_smoothed_tv_Tom_slab_107.mod"
TABLEFILE_PURE  = r"/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_pure_Tom_slab_107.mod"

WARMCOM_MODELS = [
    ("clean", TABLEFILE_CLEAN),
    ("raw",   TABLEFILE_RAW),
    ("pure",  TABLEFILE_PURE),
]

# Normalization and fit settings
E_REF = 0.15            # normalization energy (keV) relative to CLEAN
E_MAX_FIT = 1.0         # Gamma fit upper bound; lower bound is dynamic per spectrum
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

te_warm  = as_float(input_with_default("warmcom Te (keV)", "0.32"))
tau_warm = as_float(input_with_default("warmcom tau", "12.0"))
z_value  = as_float(input_with_default("redshift z", "0.0"))

# ==================== Preflight ====================
missing = [p for _, p in WARMCOM_MODELS if not os.path.exists(p)]
if missing:
    raise FileNotFoundError("Missing .mod files:\n" + "\n".join(missing))

os.makedirs(TARGET_DIR, exist_ok=True)
have_xspec = which("xspec") is not None
if not have_xspec:
    print("[ERROR] XSPEC not found in PATH.")
    sys.exit(1)

# ==================== Helpers ====================
def rm_if_exists(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception as e:
        print(f"[WARN] Cannot remove existing file {path}: {e}")

def write_xcm_warmcom(tag, mod_path, te, tau, z, out_dir):
    """
    为每个模型生成各自的 XCM，并把 QDP 输出到唯一文件：
    clean -> c.qdp, raw -> r.qdp, pure -> p.qdp
    """
    xcm = os.path.join(out_dir, f"{tag}.xcm")
    short = {"slab_clean": "c", "slab_raw": "r", "slab_pure": "p"}.get(tag, tag[:1])
    qdp_file = os.path.join(out_dir, f"{short}.qdp")
    rm_if_exists(qdp_file)

    xcm_txt = f"""
cd {out_dir}
model atable{{{mod_path}}}
{te}
{tau}
{z}
1
energies {E_MIN} {E_MAX} {NPTS} log
cpd /null
plot model
setplot command wdata {short}
setplot command exit
iplot
exit
""".lstrip()
    with open(xcm, "w") as f:
        f.write(xcm_txt)
    return xcm, qdp_file

def run_xspec(xcm):
    try:
        subprocess.run(["xspec", "-", xcm], check=True,
                       cwd=TARGET_DIR,  # 双保险：在目标目录运行
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
                    x.append(float(cols[0])); y.append(float(cols[2]))
                except Exception:
                    pass
    return np.array(x, dtype=float), np.array(y, dtype=float)

def interp_safe(x_new, x_old, y_old):
    return np.interp(x_new, x_old, y_old, left=np.nan, right=np.nan)

def first_valid_energy(E, F):
    m = np.isfinite(E) & np.isfinite(F) & (E > 0) & (F > 0)
    if not np.any(m):
        return np.nan
    return float(E[m][0])

def range_fit_loglog(E, F, Emin, Emax, min_pts=10):
    """
    log F = m * log E + b ; return Gamma = -m, together with (m, b).
    """
    mask = np.isfinite(E) & np.isfinite(F) & (E > 0) & (F > 0) & (E >= Emin) & (E <= Emax)
    E2, F2 = E[mask], F[mask]
    if E2.size < min_pts or Emin >= Emax:
        return np.nan, np.nan, np.nan
    x, y = np.log(E2), np.log(F2)
    m, b = np.polyfit(x, y, 1)
    return -m, m, b

# ==================== Run XSPEC for three models ====================
jobs = []
for name, path in WARMCOM_MODELS:
    tag = f"slab_{name}"  # slab_clean / slab_raw / slab_pure
    xcm_w, qdp_w = write_xcm_warmcom(tag, path, te_warm, tau_warm, z_value, TARGET_DIR)
    if not run_xspec(xcm_w):
        sys.exit(f"[ERROR] XSPEC failed on {name}.")
    jobs.append((name, path, tag, xcm_w, qdp_w))

# ==================== Read QDP & compute/plot ====================
# Read CLEAN first (as reference: c.qdp)
ref_job = next(j for j in jobs if j[0] == "clean")
_, _, _, _, ref_qdp = ref_job
E_ref, F_ref = read_qdp(ref_qdp)
if E_ref.size == 0:
    sys.exit(f"[ERROR] CLEAN QDP not found/readable at {ref_qdp}.")

F_ref_at_Eref = np.interp(E_REF, E_ref, F_ref, left=np.nan, right=np.nan)
if not np.isfinite(F_ref_at_Eref) or F_ref_at_Eref <= 0:
    sys.exit("[ERROR] CLEAN normalization point invalid at E_REF.")

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
pdf_name = f"cmp_warm_slab_triplet_Te{str(te_warm).replace('.','p')}_tau{str(tau_warm).replace('.','p')}_z{str(z_value).replace('.','p')}_{timestamp}.pdf"
pdf_path = os.path.join(TARGET_DIR, pdf_name)

results = []

# CLEAN (reference), fit on its native array with dynamic Emin
Emin_ref = first_valid_energy(E_ref, F_ref)
Gamma_ref, m_ref, b_ref = range_fit_loglog(E_ref, F_ref, Emin_ref, E_MAX_FIT)
print(f"[clean] fit range = [{Emin_ref:.3f}, {E_MAX_FIT:.3f}] keV | Gamma = {Gamma_ref:.4f}")

results.append({
    "label": "clean (ref)",
    "color": "C0",
    "E": E_ref,
    "F": F_ref,
    "Gamma": Gamma_ref,
    "m": m_ref,
    "b": b_ref,
    "Emin_fit": Emin_ref,
    "scaled": False,
})

# RAW & PURE: scale to CLEAN at E_REF, resample onto E_ref grid, then fit
color_map = {"raw": "C1", "pure": "C2"}
for name, _, tag, _, qdp in jobs:
    if name == "clean":
        continue
    E_m, F_m = read_qdp(qdp)
    if E_m.size == 0:
        sys.exit(f"[ERROR] Missing QDP for {name} at {qdp}.")
    F_m_at_Eref = np.interp(E_REF, E_m, F_m, left=np.nan, right=np.nan)
    if not np.isfinite(F_m_at_Eref) or F_m_at_Eref == 0:
        sys.exit(f"[ERROR] Bad normalization at E_REF for {name}.")
    scale = F_ref_at_Eref / F_m_at_Eref

    F_m_scaled_on_Eref = interp_safe(E_ref, E_m, scale * F_m)
    Emin_m = first_valid_energy(E_ref, F_m_scaled_on_Eref)
    Gamma_m, m_m, b_m = range_fit_loglog(E_ref, F_m_scaled_on_Eref, Emin_m, E_MAX_FIT)
    print(f"[{name}] scale={scale:.4g} | fit range=[{Emin_m:.3f},{E_MAX_FIT:.3f}] keV | Gamma = {Gamma_m:.4f}")

    results.append({
        "label": f"{name} × {scale:.3g}",
        "color": color_map.get(name, "C3"),
        "E": E_ref,  # already on reference grid
        "F": F_m_scaled_on_Eref,
        "Gamma": Gamma_m,
        "m": m_m,
        "b": b_m,
        "Emin_fit": Emin_m,
        "scaled": True
    })

# ====== PDF page 1: Spectra + fitted straight lines ======
param_text = (
    f"warmcom slab @ (Te={te_warm:.3f} keV, tau={tau_warm:.3f}, z={z_value}); "
    f"normalize RAW/PURE to CLEAN at {E_REF:.3f} keV\n"
    f"Dashed lines: per-spectrum log–log fits over [Emin_data, {E_MAX_FIT:.1f}] keV."
)

from matplotlib import patheffects
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

def _end_label(ax, E, F, text, color):
    m = np.isfinite(E) & np.isfinite(F) & (F > 0)
    if not np.any(m):
        return
    x, y = E[m][-1], F[m][-1]
    ax.annotate(text, xy=(x, y), xytext=(5, 0), textcoords="offset points",
                color=color, fontsize=8, va="center", ha="left",
                path_effects=[patheffects.withStroke(linewidth=2, foreground="white")])

# 统一样式：颜色/线型/marker
STYLE = {
    "clean": dict(color="k",  ls="-",  lw=2.0, marker=None,   markevery=None),
    "raw":   dict(color="C1", ls="--", lw=1.4, marker="o",    markevery=60),
    "pure":  dict(color="C2", ls="-.", lw=1.4, marker="s",    markevery=60),
}

with PdfPages(pdf_path) as pdf:
    # ===== Page 1: Spectra + fitted lines + inset zoom =====
    fig1, ax = plt.subplots()
    for r in results:
        # 找到 key
        key = "clean" if r["label"].startswith("clean") else ("raw" if r["label"].startswith("raw") else "pure")
        st = STYLE[key]
        # 主曲线
        line = ax.loglog(r["E"], r["F"], label=f"{r['label']} | Γ_fit={r['Gamma']:.3f}",
                         color=st["color"], ls=st["ls"], lw=st["lw"],
                         marker=st["marker"], markevery=st["markevery"], ms=3, alpha=0.95,
                         path_effects=[patheffects.withStroke(linewidth=1.2, foreground="white")])
        # 拟合虚线
        if np.isfinite(r["m"]) and np.isfinite(r["b"]) and np.isfinite(r["Emin_fit"]) and r["Emin_fit"] < E_MAX_FIT:
            Efit = np.geomspace(r["Emin_fit"], E_MAX_FIT, 180)
            yfit = np.exp(r["m"] * np.log(Efit) + r["b"])
            ax.loglog(Efit, yfit, ls=":", lw=1.1, color=st["color"], alpha=0.9)

        # 末端尾注
        _end_label(ax, r["E"], r["F"], key, st["color"])

    ax.axvline(E_REF, ls="--", lw=0.8, color="gray")
    ax.set_xlabel("Energy [keV]"); ax.set_ylabel("Model spectrum")
    ax.set_title(f"Three slab warmcom models (clean/raw/pure); dashed = fit over [Emin_data, {E_MAX_FIT:.1f}] keV")

    # 图例放外侧，避免压线
    leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=True)
    for t in leg.get_texts():
        t.set_path_effects([patheffects.withStroke(linewidth=2, foreground="white")])

    # 参数说明
    param_text = (
        f"Te={te_warm:.3f} keV, tau={tau_warm:.3f}, z={z_value}; "
        f"normalize RAW/PURE to CLEAN at {E_REF:.3f} keV\n"
        f"Dashed lines: per-spectrum log–log fits over [Emin_data, {E_MAX_FIT:.1f}] keV."
    )
    fig1.text(0.01, 0.98, param_text, ha="left", va="top",
              bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85), fontsize=8)

    # Inset zoom（可调整 xlim/ylim 聚焦最想看的能区）
    axins = inset_axes(ax, width="40%", height="45%", loc="lower left", borderpad=1.2)
    x1, x2 = 0.12, 0.30          # 放大能区（自己改）
    axins.set_xlim(x1, x2)
    # 自动推一个合适的 y 范围
    yvals = []
    for r in results:
        m = (r["E"] >= x1) & (r["E"] <= x2) & np.isfinite(r["F"]) & (r["F"] > 0)
        if np.any(m):
            yvals.append(r["F"][m])
    if yvals:
        ycat = np.concatenate(yvals)
        lo, hi = np.nanpercentile(ycat, [5, 95])
        axins.set_ylim(max(lo*0.8, 1e-99), hi*1.2)

    for r in results:
        key = "clean" if r["label"].startswith("clean") else ("raw" if r["label"].startswith("raw") else "pure")
        st = STYLE[key]
        axins.loglog(r["E"], r["F"], color=st["color"], ls=st["ls"], lw=st["lw"],
                     marker=st["marker"], markevery=st["markevery"], ms=2, alpha=0.95)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5", lw=0.8)

    pdf.savefig(fig1, bbox_inches="tight"); plt.close(fig1)

    # ===== Page 2: Residuals vs CLEAN (raw & pure 同页) =====
    ref_E, ref_F = results[0]["E"], results[0]["F"]
    fig2, ax2 = plt.subplots()
    for r in results[1:]:
        key = "raw" if r["label"].startswith("raw") else "pure"
        st = STYLE[key]
        ratio = r["F"] / ref_F
        ax2.semilogx(ref_E, ratio, label=f"{key} / clean", color=st["color"],
                     ls=st["ls"], lw=1.6, marker=st["marker"], markevery=70, ms=3, alpha=0.95)
    ax2.axhline(1.0, ls=":", lw=0.9, color="gray")
    ax2.set_xlabel("Energy [keV]"); ax2.set_ylabel("Ratio to clean")
    ax2.set_title("Ratio (raw/pure) relative to clean")
    ax2.legend(loc="best", fontsize=8)
    ax2.axvline(E_REF, ls="--", lw=0.8, color="gray")
    pdf.savefig(fig2, bbox_inches="tight"); plt.close(fig2)

    # ===== Page 3-4: 单独 residual + ratio（保留你原来每模型一页的细节） =====
    for r in results[1:]:
        residual = ref_F - r["F"]
        ratio    = r["F"] / ref_F
        key = "raw" if r["label"].startswith("raw") else "pure"
        st = STYLE[key]

        fig, axr = plt.subplots()
        axr.semilogx(ref_E, residual, label="Residual (clean - model×s)",
                     color=st["color"], lw=1.6)
        axr.axhline(0, ls="--", lw=0.8)
        axr.set_xlabel("Energy [keV]"); axr.set_ylabel("Residual")

        axr2 = axr.twinx()
        axr2.semilogx(ref_E, ratio, color=st["color"], ls=st["ls"], alpha=0.8,
                      label=f"Ratio ({key}/clean)")
        axr2.axhline(1.0, ls=":", lw=0.8)
        axr2.set_ylabel("Ratio")

        # 合并图例
        L1, lab1 = axr.get_legend_handles_labels()
        L2, lab2 = axr2.get_legend_handles_labels()
        axr.legend(L1 + L2, lab1 + lab2, loc="best", fontsize=8)

        info = (f"{r['label']} | Γ_fit[{r['Emin_fit']:.3f}-{E_MAX_FIT:.3f}]={r['Gamma']:.3f}\n"
                f"{param_text}")
        fig.text(0.01, 0.98, info, ha="left", va="top",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85), fontsize=8)
        axr.axvline(E_REF, ls="--", lw=0.8); axr2.axvline(E_REF, ls="--", lw=0.8)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

