import numpy as np
import os
os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
import subprocess
from matplotlib.backends.backend_pdf import PdfPages

# ========== 输入参数 ==========
te_value   = float(input("请输入 warmcom 的 te 值: "))
tau_value  = float(input("请输入 warmcom 的 tau 值: "))
z_value    = float(input("请输入 redshift z 值: "))

gamma      = float(input("请输入 nthComp 的 Gamma 值: "))
kTe_nth    = float(input("请输入 nthComp 的 kT_e 值 (keV): "))
kTbb_nth   = float(input("请输入 nthComp 的 kTbb 值 (keV):"))
# norm_nth   = float(input("请输入 nthComp 的 norm 值: "))

# ========== 输出路径：当前运行目录 ==========
current_dir = os.getcwd()
qdp_warmcom = os.path.join(current_dir, "warmcom.qdp")
qdp_nthcomp = os.path.join(current_dir, "nthcomp.qdp")
pdf_output  = os.path.join(current_dir, "compare_warmcom_nthcomp.pdf")
print(current_dir)
# ========== 生成 warmcom 脚本 ==========
xcm_warmcom = f"""
model atable{{/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-2.0_10-20_smoothed_tv_clean_Tom.mod}}
{te_value}
{tau_value}
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
1
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

# ========== 运行 XSPEC 脚本 ==========
subprocess.run(["xspec", "-", "warmcom_run.xcm"])
subprocess.run(["xspec", "-", "nthcomp_run.xcm"])

# ========== 读取 QDP 文件 ==========
def read_qdp(qdp_file):
    x, y = [], []
    with open(qdp_file, "r") as f:
        for line in f:
            if not line.startswith(("!", "@", "READ")):
                cols = line.split()
                try:
                    x.append(float(cols[0]))
                    y.append(float(cols[2]))
                except:
                    continue
    return np.array(x), np.array(y)

E_warm, F_warm = read_qdp(qdp_warmcom)
E_nth,  F_nth  = read_qdp(qdp_nthcomp)



# —— 计算 nthComp 的 tau（以 nthComp 为标准的近似关系）——
me_keV = 511.0
theta_e = kTe_nth / me_keV
term = (gamma + 0.5)**2 - 2.25
if term <= 0 or theta_e <= 0:
    raise RuntimeError("参数不物理：需满足 (Gamma+0.5)^2 > 9/4 且 kT_e > 0。")
tau_nth = (1.0 / (theta_e * term))**0.5 - 1.0/3.0
print(f"[nthComp] 由 Γ={gamma:.3f}, kT_e={kTe_nth:.3f} keV 推得 tau ≈ {tau_nth:.3f}")

# ========== 归一化（比如归一化到1 keV 附近） ==========
# ===== 选择最小二乘对齐的能段 =====
E1, E2 = 0.3, 2.0   # 你可改，比如 (0.5, 10.0)

# —— 若能量轴不同，把 nthComp 插值到 warmcom 的能量点 ——
if not np.array_equal(E_warm, E_nth):
    F_nth_on_warm = np.interp(E_warm, E_nth, F_nth, left=np.nan, right=np.nan)
else:
    F_nth_on_warm = F_nth.copy()

# —— 选取对齐能段上的点，并去掉 NaN/无效点 ——
m_band = (E_warm >= E1) & (E_warm <= E2) & np.isfinite(F_warm) & np.isfinite(F_nth_on_warm)
if not np.any(m_band):
    raise RuntimeError(f"在对齐能段 [{E1},{E2}] keV 没有有效数据点；请调整 E1/E2。")

Fw = F_warm[m_band]
Fn = F_nth_on_warm[m_band]

# ===== 以 nthComp 为标准的最小二乘尺度因子：min || s*Fw - Fn ||^2 =====
den = np.sum(Fw*Fw)
if den == 0 or not np.isfinite(den):
    raise RuntimeError("warmcom 在所选能段的谱为零或无效，无法归一化。")
s = np.sum(Fw*Fn) / den

# —— 计算全能段下的缩放后 warmcom、残差与比例 ——
Fw_scaled = s * F_warm
residual  = F_nth_on_warm - Fw_scaled          # 绝对残差（同单位）
ratio     = Fw_scaled / F_nth_on_warm          # scaled_warm / nthComp
m_valid   = np.isfinite(E_warm) & np.isfinite(Fw_scaled) & np.isfinite(F_nth_on_warm) & (F_nth_on_warm > 0)

# ===== 保存到同一个 PDF（两页）=====
with PdfPages(pdf_output) as pdf:
    # 第1页：对比（nthComp vs s*warmcom）
    fig1 = plt.figure()
    plt.loglog(E_warm[m_valid], F_nth_on_warm[m_valid],  label="nthComp (reference)")
    plt.loglog(E_warm[m_valid], Fw_scaled[m_valid],      label=f"warmcom × s (s={s:.4g})")
    plt.xlabel("Energy [keV]")
    plt.ylabel("Spectrum")
    plt.legend()
    pdf.savefig(fig1); plt.close(fig1)

    # 第2页：残差与/或比例（你可二选一；下面我两条都画在同页的双轴）
    fig2, ax1 = plt.subplots()
    ax1.semilogx(E_warm[m_valid], residual[m_valid], label="Residual (nth - s*warm)", lw=1.2)
    ax1.axhline(0.0, ls="--", lw=0.8, color="gray")
    ax1.set_xlabel("Energy [keV]")
    ax1.set_ylabel("Residual")

    ax2 = ax1.twinx()
    ax2.semilogx(E_warm[m_valid], ratio[m_valid], label="Ratio (s*warm / nth)", lw=1.0, alpha=0.75)
    ax2.axhline(1.0, ls=":", lw=0.8, color="gray")
    ax2.set_ylabel("Ratio")

    # 合并图例
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines+lines2, labels+labels2, loc="best")

    pdf.savefig(fig2); plt.close(fig2)

print(f"[OK] 最小二乘对齐完成：s = {s:.6g}（能段 [{E1},{E2}] keV）")
print(f"[OK] 对比/残差已保存：{pdf_output}")