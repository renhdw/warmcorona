#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
import subprocess
from matplotlib.backends.backend_pdf import PdfPages

def input_with_default(prompt, default_str):
    s = input(f"{prompt} [默认 {default_str}]: ").strip()
    return s if s else default_str

# ========== 输入：nth 的 kTe、Gamma，再给 z（kTbb 可选）==========
kTe_nth   = float(input("请输入 nthComp 的 kT_e (keV): ").strip())
gamma     = float(input("请输入 nthComp 的 Gamma: ").strip())
z_value   = float(input("请输入 redshift z 值: ").strip())
kTbb_nth  = float(input_with_default("请输入 nthComp 的 kTbb (keV)", "0.003"))

# ========== 由 (Gamma, kTe) 计算 tau (nth 的近似关系) ==========
me_keV  = 511.0
theta_e = kTe_nth / me_keV
term_g  = (gamma + 0.5)**2 - 2.25  # 2.25 = 9/4
if theta_e <= 0 or term_g <= 0:
    raise RuntimeError("参数不物理：需满足 kT_e>0 且 (Gamma+0.5)^2 > 9/4。")

tau_nth = (1.0 / (theta_e * term_g))**0.5 - 1.0/3.0
print(f"[nthComp] 由 Γ={gamma:.3f}, kT_e={kTe_nth:.3f} keV 推得 tau ≈ {tau_nth:.4f}")

# warmcom 使用同一 te 与该 tau
te_warm   = kTe_nth
tau_warm  = tau_nth

# ========== 路径 ==========
current_dir = os.getcwd()
qdp_warmcom = os.path.join(current_dir, "warmcom.qdp")
qdp_nthcomp = os.path.join(current_dir, "nthcomp.qdp")
pdf_output  = os.path.join(current_dir, "compare_warmcom_nthcomp.pdf")
print("[路径]", current_dir)

# ========== 生成 warmcom 脚本（按你的模型路径改）==========
xcm_warmcom = f"""
model atable{{/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-1.0_10-20_smoothed_tv_clean_Tom_slab.mod}}
{te_warm}
{tau_warm}
{z_value}
1
energies 0.1 10 1000 log
cpd /null
plot model
setplot command wdata {qdp_warmcom}
setplot command exit
iplot
exit
"""
with open("warmcom_run.xcm", "w") as f:
    f.write(xcm_warmcom)

# ========== 生成 nthComp 脚本 ==========
xcm_nthcomp = f"""
model nthComp
{gamma}
{max(kTe_nth, 0.01)}
{kTbb_nth}
0
{z_value}
1
newpar 2 {kTe_nth} 1e-3 0.1 0.1 10 1000
cpd /null
energies 0.1 10 1000 log
plot model
setplot command wdata {qdp_nthcomp}
setplot command exit
iplot
exit
"""
with open("nthcomp_run.xcm", "w") as f:
    f.write(xcm_nthcomp)

# ========== 运行 XSPEC ==========
subprocess.run(["xspec", "-", "warmcom_run.xcm"])
subprocess.run(["xspec", "-", "nthcomp_run.xcm"])

# ========== 读取 QDP ==========
def read_qdp(qdp_file):
    x, y = [], []
    with open(qdp_file, "r") as f:
        for line in f:
            if not line.startswith(("!", "@", "READ")):
                cols = line.split()
                if len(cols) >= 3:
                    try:
                        x.append(float(cols[0]))   # E
                        y.append(float(cols[2]))   # model
                    except:
                        pass
    return np.array(x, dtype=float), np.array(y, dtype=float)

E_warm, F_warm = read_qdp(qdp_warmcom)
E_nth,  F_nth  = read_qdp(qdp_nthcomp)

# ========== 在 E_ref 归一化（单点）=========
E_ref = 0.15  # keV
F_w_ref = np.interp(E_ref, E_warm, F_warm, left=np.nan, right=np.nan)
F_n_ref = np.interp(E_ref, E_nth,  F_nth,  left=np.nan, right=np.nan)
if not np.isfinite(F_w_ref) or not np.isfinite(F_n_ref) or F_w_ref == 0:
    raise RuntimeError("无法在 0.15 keV 处进行归一化：请检查能量网格或模型输出。")
s = F_n_ref / F_w_ref

# 统一能量轴用于对比与残差
if not np.array_equal(E_warm, E_nth):
    F_nth_on_warm = np.interp(E_warm, E_nth, F_nth, left=np.nan, right=np.nan)
else:
    F_nth_on_warm = F_nth.copy()

Fw_scaled = s * F_warm
residual  = F_nth_on_warm - Fw_scaled
ratio     = Fw_scaled / F_nth_on_warm
m_valid   = np.isfinite(E_warm) & np.isfinite(Fw_scaled) & np.isfinite(F_nth_on_warm) & (F_nth_on_warm > 0)

# ========== 新增：对 warm 谱做局部 log–log 线性拟合，得到 Γ_warm，并生成拟合曲线 ==========
def local_gamma(E, F, E0, factor=1.5, min_pts=7):
    mask = np.isfinite(E) & np.isfinite(F) & (E > 0) & (F > 0)
    E, F = E[mask], F[mask]
    if len(E) < 2:
        return np.nan, None, None
    sel = (E >= E0 / factor) & (E <= E0 * factor)
    if np.sum(sel) >= min_pts:
        idx = np.where(sel)[0]
    else:
        k = min(len(E), min_pts)
        idx = np.argsort(np.abs(E - E0))[:k]
    x = np.log(E[idx])
    y = np.log(F[idx])
    if len(x) < 2:
        return np.nan, None, None
    m, b = np.polyfit(x, y, 1)   # y = m x + b
    return -m, m, b             # 返回 Γ, 斜率, 截距

Gamma_warm_local, m_fit, b_fit = local_gamma(E_warm, F_warm, E_ref)
print(f"[warmcom] 局部谱指数 Γ_warm(E≈{E_ref} keV) ≈ {Gamma_warm_local:.4f}")

# 生成拟合曲线（在 log-log 下是直线）
if m_fit is not None:
    F_warm_fit = np.exp(m_fit * np.log(E_warm) + b_fit)
else:
    F_warm_fit = None

# ========== 生成 PDF ==========
param_text = (
    "warmcom vs nthComp (te,tau from nth inputs)\n"
    f"z={z_value:.4g} | nth: kTe={kTe_nth:.4f} keV, Gamma={gamma:.4f} -> tau={tau_nth:.4f}\n"
    f"warmcom: te={te_warm:.4f} keV, tau={tau_warm:.4f} | Γ_warm(local @ {E_ref:.3f} keV)={Gamma_warm_local:.4f}\n"
    f"normalize at {E_ref:.2f} keV | scale s={s:.6g}"
)

with PdfPages(pdf_output) as pdf:
    # 页1：模型对比
    fig1 = plt.figure()
    plt.loglog(E_warm[m_valid], F_nth_on_warm[m_valid],  label="nthComp (input Γ,kTe)")
    plt.loglog(E_warm[m_valid], Fw_scaled[m_valid],      label="warmcom × s (E_ref normalized)")
    if F_warm_fit is not None:
        plt.loglog(E_warm[m_valid], s*F_warm_fit[m_valid], '--', color='red',
                   label="warmcom local fit")
    plt.xlabel("Energy [keV]")
    plt.ylabel("Spectrum")
    plt.legend(loc="best")
    fig1.text(0.01, 0.98, param_text, ha="left", va="top",
              bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.8), fontsize=9)
    plt.axvline(E_ref, ls="--", lw=0.8)
    fig1.text(0.99, 0.02, f"Γ_warm(local) ≈ {Gamma_warm_local:.3f}",
              ha="right", va="bottom",
              bbox=dict(boxstyle="round,pad=0.35", fc="white", alpha=0.85), fontsize=9)
    pdf.savefig(fig1); plt.close(fig1)

    # 页2：残差与比例
    fig2, ax1 = plt.subplots()
    ax1.semilogx(E_warm[m_valid], residual[m_valid], label="Residual (nth - s*warm)", lw=1.2)
    ax1.axhline(0.0, ls="--", lw=0.8)
    ax1.set_xlabel("Energy [keV]")
    ax1.set_ylabel("Residual")
    ax2 = ax1.twinx()
    ax2.semilogx(E_warm[m_valid], ratio[m_valid], label="Ratio (s*warm / nth)", lw=1.0, alpha=0.75)
    ax2.axhline(1.0, ls=":", lw=0.8)
    ax2.set_ylabel("Ratio")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines+lines2, labels+labels2, loc="best")
    fig2.text(0.01, 0.98, param_text, ha="left", va="top",
              bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.8), fontsize=9)
    ax1.axvline(E_ref, ls="--", lw=0.8)
    ax2.axvline(E_ref, ls="--", lw=0.8)
    pdf.savefig(fig2); plt.close(fig2)

print(f"[OK] 已在 {E_ref} keV 归一化，s = {s:.6g}")
print(f"[OK] Γ_warm(local @ {E_ref} keV) ≈ {Gamma_warm_local:.6f}")
print(f"[OK] 对比/残差已保存：{pdf_output}")
