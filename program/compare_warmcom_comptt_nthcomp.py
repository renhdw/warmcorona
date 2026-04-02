#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import uuid
import shutil
import subprocess
import numpy as np
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# =========================================================
# 用户参数：只改这里
# =========================================================

# ---------- 物理参数 ----------
KTE   = 0.55          # keV
GAMMA = 2.8
KTBB  = 0.003        # keV = 3 eV
Z     = 0.0
NORM  = 1.0

# ---------- 输出目录 ----------
OUT_DIR = os.path.expanduser("~/data/monk/plot/warmcorona/test/2026.03.12")

# ---------- warmcom sphere ----------
SPHERE_MODEL_DIR  = os.path.expanduser("~/data/monk/plot/warmcorona/model/warmcom/warmcom_sphere_log/warmcom_sphere")
SPHERE_PKG        = "warmcom_sphere"
SPHERE_MODEL_NAME = "warmcomsphere"
SPHERE_BASEDIR    = os.path.expanduser("~/data/monk/plot/warmcorona/test/test_smooth_10_7/_exports")

# ---------- warmcom slab ----------
SLAB_MODEL_DIR  = os.path.expanduser("~/data/monk/plot/warmcorona/model/warmcom/warmcom_slab_log/warmcom_slab")
SLAB_PKG        = "warmcom_slab"
SLAB_MODEL_NAME = "warmcomslab"
SLAB_BASEDIR    = os.path.expanduser("~/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports")

# ---------- compTT geometry switch ----------
# compTT 手册：
#   par5 <= 1  -> disk
#   par5 >  1  -> sphere
COMPTT_SWITCH_SPHERE = 2.0
COMPTT_SWITCH_SLAB   = 1.0

# ---------- 能量网格 / 绘图 ----------
E_MIN     = 0.01
E_MAX     = 20.0
NBINS     = 1000
GRID_MODE = "log"

# 归一化参考能量
E_REF = 0.1

# 第三页软能段对比上限
E_SOFT_MAX = 2.0

# y 轴下限
FLUX_YMIN = 1e-6


# =========================================================
# 小工具
# =========================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def must_exist_dir(path, name):
    if not os.path.isdir(path):
        raise RuntimeError(f"{name} 不存在：{path}")

def rm_if_exists(path):
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception:
        pass

def read_qdp_model(path):
    xs, ys = [], []
    if not os.path.exists(path):
        raise RuntimeError(f"QDP 不存在：{path}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if (not s) or s[0] in ("!", "@"):
                continue
            parts = s.split()
            try:
                if len(parts) >= 3:
                    x = float(parts[0])
                    y = float(parts[2])
                elif len(parts) >= 2:
                    x = float(parts[0])
                    y = float(parts[1])
                else:
                    continue
            except Exception:
                continue
            xs.append(x)
            ys.append(y)

    if not xs:
        raise RuntimeError(f"QDP 为空：{path}")

    return np.array(xs, float), np.array(ys, float)

def interp_safe(x_new, x_old, y_old):
    return np.interp(x_new, x_old, y_old, left=np.nan, right=np.nan)

def compute_abs_err(y_model, y_ref):
    y_model = np.asarray(y_model, float)
    y_ref   = np.asarray(y_ref, float)
    ok = np.isfinite(y_model) & np.isfinite(y_ref) & (y_model > 0) & (y_ref > 0)
    out = np.full_like(y_model, np.nan, dtype=float)
    out[ok] = np.abs(1.0 - y_model[ok] / y_ref[ok])
    return out

def compute_ratio(y_model, y_ref):
    y_model = np.asarray(y_model, float)
    y_ref   = np.asarray(y_ref, float)
    ok = np.isfinite(y_model) & np.isfinite(y_ref) & (y_model > 0) & (y_ref > 0)
    out = np.full_like(y_model, np.nan, dtype=float)
    out[ok] = y_model[ok] / y_ref[ok]
    return out


# =========================================================
# 由 nthComp 的 (Gamma, kTe) 推 tau
# sphere: Zdziarski et al. (1996) 形式
# slab  : Titarchuk (1994) 近似公式
# =========================================================

def gamma_kte_to_tau_sphere(gamma, kte_keV):
    me_keV = 511.0
    theta_e = kte_keV / me_keV
    term = (gamma + 0.5)**2 - 2.25   # 9/4
    if theta_e <= 0 or term <= 0:
        raise RuntimeError(
            f"Non-physical sphere parameters: kTe={kte_keV}, Gamma={gamma}"
        )
    tau = np.sqrt(2.25 + 3.0 / (theta_e * term)) - 1.5
    return float(tau)

def gamma_kte_to_tau_slab(gamma, kte_keV):
    me_keV = 511.0
    theta_e = kte_keV / me_keV
    term = (gamma + 0.5)**2 - 2.25   # 9/4
    if theta_e <= 0 or term <= 0:
        raise RuntimeError(
            f"Non-physical slab parameters: kTe={kte_keV}, Gamma={gamma}"
        )
    tau = np.sqrt(1.0 / (theta_e * term)) - 1.0 / 3.0
    return float(tau)


# =========================================================
# XSPEC 跑模型并导出 QDP
# =========================================================

def run_xspec_with_cmds(cmds, out_dir, env_extra=None):
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env.setdefault("PGPLOT_DEV", "/null")
    if env_extra:
        env.update(env_extra)

    proc = subprocess.run(
        ["xspec"],
        input=cmds,
        text=True,
        cwd=out_dir,
        env=env,
        capture_output=True
    )

    if proc.returncode != 0:
        print("\n[XSPEC STDOUT tail]\n", proc.stdout[-4000:])
        print("\n[XSPEC STDERR tail]\n", proc.stderr[-4000:])
        raise RuntimeError(f"XSPEC failed, returncode={proc.returncode}")

def write_qdp_warmcom(model_dir, pkg, model_name, basedir,
                      te, tau, z, norm, out_qdp, out_xcm):
    ensure_dir(os.path.dirname(out_qdp))
    tmp_short = f"tmp_{uuid.uuid4().hex[:8]}.qdp"
    tmp_path  = os.path.join(os.path.dirname(out_qdp), tmp_short)
    rm_if_exists(tmp_path)
    rm_if_exists(out_qdp)

    grid_kw = "log" if GRID_MODE.lower().startswith("log") else "lin"

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
energies {E_MIN} {E_MAX} {NBINS} {grid_kw}
plot model
iplot
wdata {tmp_short}
quit
exit
"""

    with open(out_xcm, "w", encoding="utf-8") as f:
        f.write(cmds)

    run_xspec_with_cmds(cmds, os.path.dirname(out_qdp), env_extra={"WARMCOM_BASEDIR": basedir})

    if (not os.path.exists(tmp_path)) or os.path.getsize(tmp_path) == 0:
        raise RuntimeError(f"warmcom 未生成 QDP：{tmp_path}")

    shutil.move(tmp_path, out_qdp)

def write_qdp_nthcomp(gamma, kte, ktbb, z, out_qdp, out_xcm):
    ensure_dir(os.path.dirname(out_qdp))
    tmp_short = f"tmp_{uuid.uuid4().hex[:8]}.qdp"
    tmp_path  = os.path.join(os.path.dirname(out_qdp), tmp_short)
    rm_if_exists(tmp_path)
    rm_if_exists(out_qdp)

    grid_kw = "log" if GRID_MODE.lower().startswith("log") else "lin"
    kte_plot = max(float(kte), 0.01)

    cmds = f"""query yes
cpd /null
model nthComp
{gamma}
{kte_plot}
{ktbb}
0
{z}
1
newpar 2 {kte} 1e-3 0.1 0.1 10 1000
setplot rebin 1 1
setplot energy
setplot area off
energies {E_MIN} {E_MAX} {NBINS} {grid_kw}
plot model
iplot
wdata {tmp_short}
quit
exit
"""

    with open(out_xcm, "w", encoding="utf-8") as f:
        f.write(cmds)

    run_xspec_with_cmds(cmds, os.path.dirname(out_qdp))

    if (not os.path.exists(tmp_path)) or os.path.getsize(tmp_path) == 0:
        raise RuntimeError(f"nthComp 未生成 QDP：{tmp_path}")

    shutil.move(tmp_path, out_qdp)

def write_qdp_comptt(z, t0, kte, tau, geom_switch, norm, out_qdp, out_xcm):
    ensure_dir(os.path.dirname(out_qdp))
    tmp_short = f"tmp_{uuid.uuid4().hex[:8]}.qdp"
    tmp_path  = os.path.join(os.path.dirname(out_qdp), tmp_short)
    rm_if_exists(tmp_path)
    rm_if_exists(out_qdp)

    grid_kw = "log" if GRID_MODE.lower().startswith("log") else "lin"

    cmds = f"""query yes
cpd /null
model compTT
{z}
{t0}
{kte}
{tau}
{geom_switch}
{norm}
freeze 5
newpar 3 {kte} 1e-3 0.1 0.1 10 1000
setplot rebin 1 1
setplot energy
setplot area off
energies {E_MIN} {E_MAX} {NBINS} {grid_kw}
plot model
iplot
wdata {tmp_short}
quit
exit
"""

    with open(out_xcm, "w", encoding="utf-8") as f:
        f.write(cmds)

    run_xspec_with_cmds(cmds, os.path.dirname(out_qdp))

    if (not os.path.exists(tmp_path)) or os.path.getsize(tmp_path) == 0:
        raise RuntimeError(f"compTT 未生成 QDP：{tmp_path}")

    shutil.move(tmp_path, out_qdp)


# =========================================================
# 统一到 nthComp 网格，并在 E_REF 归一化
# =========================================================

def prepare_page_results(E_ref, F_ref, others):
    """
    E_ref, F_ref: nthComp
    others: list of dict(label, E, F, color)
    """
    results = []

    results.append({
        "label": "nthComp",
        "base_label": "nthComp",
        "E": E_ref,
        "F": F_ref,
        "color": "C0",
        "scale": 1.0,
        "is_ref": True,
    })

    f_ref_at_eref = np.interp(E_REF, E_ref, F_ref, left=np.nan, right=np.nan)
    if not np.isfinite(f_ref_at_eref) or f_ref_at_eref <= 0:
        raise RuntimeError(f"参考模型 nthComp 在 E_REF={E_REF} keV 处无有效值。")

    for item in others:
        E0 = item["E"]
        F0 = item["F"]
        f0_at_eref = np.interp(E_REF, E0, F0, left=np.nan, right=np.nan)
        if not np.isfinite(f0_at_eref) or f0_at_eref <= 0:
            raise RuntimeError(f"{item['label']} 在 E_REF={E_REF} keV 处无有效值。")

        scale = f_ref_at_eref / f0_at_eref
        F_on_ref = interp_safe(E_ref, E0, scale * F0)

        results.append({
            "label": f"{item['label']} × {scale:.3g}",
            "base_label": item["label"],
            "E": E_ref,
            "F": F_on_ref,
            "color": item["color"],
            "scale": scale,
            "is_ref": False,
        })

    return results


# =========================================================
# 画图
# =========================================================

def plot_page(pdf, results, title, subtitle):
    ref = results[0]
    E = ref["E"]
    Fref = ref["F"]

    fig, axes = plt.subplots(
        3, 1, figsize=(9.0, 10.0),
        gridspec_kw={"height_ratios": [2.4, 1.2, 1.2]}
    )
    ax1, ax2, ax3 = axes

    # 1) 光谱
    for r in results:
        lw = 1.8 if r["is_ref"] else 1.3
        ax1.loglog(
            r["E"], r["F"],
            lw=lw, color=r["color"],
            label=r["label"]
        )

    ax1.axvline(E_REF, ls="--", lw=0.8, color="gray")
    ax1.set_xlim(E_MIN, E_MAX)
    ax1.set_ylim(bottom=FLUX_YMIN)
    ax1.set_ylabel("Flux")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8)
    ax1.set_title(title)

    # 2) abs err
    for r in results[1:]:
        err = compute_abs_err(r["F"], Fref)
        ax2.semilogx(E, err, lw=1.2, color=r["color"], label=r["base_label"])
    ax2.axhline(0.0, ls="--", lw=0.8, color="gray")
    ax2.axvline(E_REF, ls="--", lw=0.8, color="gray")
    ax2.set_xlim(E_MIN, E_MAX)
    ax2.set_xlabel("Energy [keV]")
    ax2.set_ylabel(r"$\left|1-\frac{\mathrm{MODEL}}{\mathrm{REF}}\right|$")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=8)

    # 3) ratio
    for r in results[1:]:
        ratio = compute_ratio(r["F"], Fref)
        ax3.semilogx(E, ratio, lw=1.2, color=r["color"], label=r["base_label"])
    ax3.axhline(1.0, ls="--", lw=0.8, color="gray")
    ax3.axvline(E_REF, ls="--", lw=0.8, color="gray")
    ax3.set_xlim(E_MIN, E_MAX)
    ax3.set_xlabel("Energy [keV]")
    ax3.set_ylabel("MODEL / REF")
    ax3.grid(True, which="both", alpha=0.25)
    ax3.legend(fontsize=8)

    fig.text(
        0.01, 0.985, subtitle,
        ha="left", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.90)
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_soft_page(pdf, lines_5, title, subtitle):
    fig, ax = plt.subplots(figsize=(9.0, 7.0))

    for r in lines_5:
        ax.loglog(r["E"], r["F"], lw=1.5, color=r["color"], label=r["label"])

    ax.axvline(E_REF, ls="--", lw=0.8, color="gray")
    ax.set_xlim(E_MIN, E_SOFT_MAX)
    ax.set_ylim(bottom=FLUX_YMIN)
    ax.set_xlabel("Energy [keV]")
    ax.set_ylabel("Flux")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)

    fig.text(
        0.01, 0.985, subtitle,
        ha="left", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.90)
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# =========================================================
# 主程序
# =========================================================

def main():
    ensure_dir(OUT_DIR)

    if shutil.which("xspec") is None:
        raise RuntimeError("找不到 xspec，请先确保 XSPEC 在 PATH 中。")

    must_exist_dir(SPHERE_MODEL_DIR, "SPHERE_MODEL_DIR")
    must_exist_dir(SPHERE_BASEDIR,   "SPHERE_BASEDIR")
    must_exist_dir(SLAB_MODEL_DIR,   "SLAB_MODEL_DIR")
    must_exist_dir(SLAB_BASEDIR,     "SLAB_BASEDIR")

    tau_sphere = gamma_kte_to_tau_sphere(GAMMA, KTE)
    tau_slab   = gamma_kte_to_tau_slab(GAMMA, KTE)

    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    tag = f"cmp_warmcom_comptt_nthcomp_kte{KTE:g}_g{GAMMA:g}_z{Z:g}_{ts}"

    out_pdf = os.path.join(OUT_DIR, f"{tag}.pdf")
    out_txt = os.path.join(OUT_DIR, f"{tag}_summary.txt")

    print("\n================ INPUT ================")
    print(f"KTE         = {KTE} keV")
    print(f"GAMMA       = {GAMMA}")
    print(f"KTBB        = {KTBB} keV")
    print(f"Z           = {Z}")
    print(f"NORM        = {NORM}")
    print(f"tau_sphere  = {tau_sphere:.8g}   [Zdziarski et al. 1996]")
    print(f"tau_slab    = {tau_slab:.8g}   [Titarchuk 1994 approx.]")
    print(f"E range     = {E_MIN} - {E_MAX} keV")
    print(f"soft range  = {E_MIN} - {E_SOFT_MAX} keV")
    print(f"OUT         = {OUT_DIR}")
    print("=======================================\n")

    # ---------------- Page 1: sphere ----------------
    qdp_nth_s = os.path.join(OUT_DIR, f"{tag}__p01_nthcomp.qdp")
    xcm_nth_s = os.path.join(OUT_DIR, f"{tag}__p01_nthcomp.xcm")

    qdp_wc_s  = os.path.join(OUT_DIR, f"{tag}__p01_warmcom_sphere.qdp")
    xcm_wc_s  = os.path.join(OUT_DIR, f"{tag}__p01_warmcom_sphere.xcm")

    qdp_ct_s  = os.path.join(OUT_DIR, f"{tag}__p01_comptt_sphere.qdp")
    xcm_ct_s  = os.path.join(OUT_DIR, f"{tag}__p01_comptt_sphere.xcm")

    print("[Page 1] sphere ...")
    write_qdp_nthcomp(GAMMA, KTE, KTBB, Z, qdp_nth_s, xcm_nth_s)
    write_qdp_warmcom(SPHERE_MODEL_DIR, SPHERE_PKG, SPHERE_MODEL_NAME, SPHERE_BASEDIR,
                      KTE, tau_sphere, Z, NORM, qdp_wc_s, xcm_wc_s)
    write_qdp_comptt(Z, KTBB, KTE, tau_sphere, COMPTT_SWITCH_SPHERE, NORM, qdp_ct_s, xcm_ct_s)

    En_s, Fn_s = read_qdp_model(qdp_nth_s)
    Ew_s, Fw_s = read_qdp_model(qdp_wc_s)
    Ec_s, Fc_s = read_qdp_model(qdp_ct_s)

    results_sphere = prepare_page_results(
        En_s, Fn_s,
        [
            {"label": "warmcom_sphere", "E": Ew_s, "F": Fw_s, "color": "C1"},
            {"label": "compTT_sphere",  "E": Ec_s, "F": Fc_s, "color": "C2"},
        ]
    )

    # ---------------- Page 2: slab ----------------
    qdp_nth_l = os.path.join(OUT_DIR, f"{tag}__p02_nthcomp.qdp")
    xcm_nth_l = os.path.join(OUT_DIR, f"{tag}__p02_nthcomp.xcm")

    qdp_wc_l  = os.path.join(OUT_DIR, f"{tag}__p02_warmcom_slab.qdp")
    xcm_wc_l  = os.path.join(OUT_DIR, f"{tag}__p02_warmcom_slab.xcm")

    qdp_ct_l  = os.path.join(OUT_DIR, f"{tag}__p02_comptt_slab.qdp")
    xcm_ct_l  = os.path.join(OUT_DIR, f"{tag}__p02_comptt_slab.xcm")

    print("[Page 2] slab ...")
    write_qdp_nthcomp(GAMMA, KTE, KTBB, Z, qdp_nth_l, xcm_nth_l)
    write_qdp_warmcom(SLAB_MODEL_DIR, SLAB_PKG, SLAB_MODEL_NAME, SLAB_BASEDIR,
                      KTE, tau_slab, Z, NORM, qdp_wc_l, xcm_wc_l)
    write_qdp_comptt(Z, KTBB, KTE, tau_slab, COMPTT_SWITCH_SLAB, NORM, qdp_ct_l, xcm_ct_l)

    En_l, Fn_l = read_qdp_model(qdp_nth_l)
    Ew_l, Fw_l = read_qdp_model(qdp_wc_l)
    Ec_l, Fc_l = read_qdp_model(qdp_ct_l)

    results_slab = prepare_page_results(
        En_l, Fn_l,
        [
            {"label": "warmcom_slab", "E": Ew_l, "F": Fw_l, "color": "C3"},
            {"label": "compTT_disk",  "E": Ec_l, "F": Fc_l, "color": "C4"},
        ]
    )

    # ---------------- Page 3: soft-band 5 lines ----------------
    soft_lines = [
        {
            "label": results_sphere[0]["label"],   # nthComp
            "E": results_sphere[0]["E"],
            "F": results_sphere[0]["F"],
            "color": results_sphere[0]["color"],
        },
        {
            "label": results_sphere[1]["label"],   # warmcom_sphere
            "E": results_sphere[1]["E"],
            "F": results_sphere[1]["F"],
            "color": results_sphere[1]["color"],
        },
        {
            "label": results_sphere[2]["label"],   # compTT_sphere
            "E": results_sphere[2]["E"],
            "F": results_sphere[2]["F"],
            "color": results_sphere[2]["color"],
        },
        {
            "label": results_slab[1]["label"],     # warmcom_slab
            "E": results_slab[1]["E"],
            "F": results_slab[1]["F"],
            "color": results_slab[1]["color"],
        },
        {
            "label": results_slab[2]["label"],     # compTT_disk
            "E": results_slab[2]["E"],
            "F": results_slab[2]["F"],
            "color": results_slab[2]["color"],
        },
    ]

    # ---------------- PDF ----------------
    with PdfPages(out_pdf) as pdf:
        subtitle1 = (
            f"Page 1 / sphere\n"
            f"kTe={KTE:g} keV; Gamma={GAMMA:g}; kTbb={KTBB:g} keV; z={Z:g}; norm={NORM:g}\n"
            f"tau_sphere={tau_sphere:.8g} [Zdziarski et al. 1996]; E_REF={E_REF:g} keV\n"
            f"energy shown = {E_MIN:g} - {E_MAX:g} keV\n"
            f"warmcom sphere: pkg={SPHERE_PKG}, model={SPHERE_MODEL_NAME}\n"
            f"compTT sphere par5={COMPTT_SWITCH_SPHERE:g}; REF=nthComp"
        )
        plot_page(pdf, results_sphere,
                  "warmcom_sphere vs nthComp vs compTT(sphere)",
                  subtitle1)

        subtitle2 = (
            f"Page 2 / slab\n"
            f"kTe={KTE:g} keV; Gamma={GAMMA:g}; kTbb={KTBB:g} keV; z={Z:g}; norm={NORM:g}\n"
            f"tau_slab={tau_slab:.8g} [Titarchuk 1994 approx.]; E_REF={E_REF:g} keV\n"
            f"energy shown = {E_MIN:g} - {E_MAX:g} keV\n"
            f"warmcom slab: pkg={SLAB_PKG}, model={SLAB_MODEL_NAME}\n"
            f"compTT disk par5={COMPTT_SWITCH_SLAB:g}; REF=nthComp"
        )
        plot_page(pdf, results_slab,
                  "warmcom_slab vs nthComp vs compTT(disk)",
                  subtitle2)

        subtitle3 = (
            f"Page 3 / soft-band summary (5 lines)\n"
            f"kTe={KTE:g} keV; Gamma={GAMMA:g}; kTbb={KTBB:g} keV; z={Z:g}; norm={NORM:g}\n"
            f"tau_sphere={tau_sphere:.8g}; tau_slab={tau_slab:.8g}; E_REF={E_REF:g} keV\n"
            f"soft energy shown = {E_MIN:g} - {E_SOFT_MAX:g} keV\n"
            f"Lines: nthComp + warmcom_sphere + compTT_sphere + warmcom_slab + compTT_disk"
        )
        plot_soft_page(pdf, soft_lines,
                       "Soft-band comparison: 5 models",
                       subtitle3)

    # ---------------- summary ----------------
    lines = []
    lines.append("# summary")
    lines.append(f"# Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("[input]")
    lines.append(f"KTE   = {KTE}")
    lines.append(f"GAMMA = {GAMMA}")
    lines.append(f"KTBB  = {KTBB}")
    lines.append(f"Z     = {Z}")
    lines.append(f"NORM  = {NORM}")
    lines.append(f"TAU_SPHERE = {tau_sphere}")
    lines.append(f"TAU_SLAB   = {tau_slab}")
    lines.append(f"E_MIN = {E_MIN}")
    lines.append(f"E_MAX = {E_MAX}")
    lines.append(f"E_SOFT_MAX = {E_SOFT_MAX}")
    lines.append("")
    lines.append("[sphere]")
    for r in results_sphere:
        lines.append(f"{r['base_label']}: scale={r['scale']:.16g}")
    lines.append("")
    lines.append("[slab]")
    for r in results_slab:
        lines.append(f"{r['base_label']}: scale={r['scale']:.16g}")
    lines.append("")
    lines.append("[outputs]")
    lines.append(f"PDF = {out_pdf}")
    lines.append(f"TXT = {out_txt}")

    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print("\n✅ DONE")
    print("PDF =", out_pdf)
    print("TXT =", out_txt)


if __name__ == "__main__":
    main()