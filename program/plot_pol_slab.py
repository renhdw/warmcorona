#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


# =========================
# 路径设置
# =========================
DATA_DIR = Path(
    "/home/hdw/data/monk/plot/warmcorona/test/polarization/data/10_7/te_0.300_tau_15.000/calspec"
)

OUT_DIR = Path(
    "~/data/monk/plot/warmcorona/test/2026.03.18/pol_107/te_0.300_tau_15.000"
).expanduser()

OUT_DIR.mkdir(parents=True, exist_ok=True)


# =========================
# 基本参数
# =========================
NBIN = 60                  # log 能量 bin 数
FLUX_FRAC = 1e-6           # 屏蔽过低通量 bin
OUT_PREFIX = "slab_binned"


# =========================
# 二进制读取函数
# =========================
def load_dat_file(file_path: Path) -> np.ndarray:
    with open(file_path, "rb") as f:
        data = f.read()

    if len(data) % 8 != 0:
        raise ValueError(f"{file_path}: 字节长度不是 8 的整数倍，无法按 double 解析。")

    arr = np.frombuffer(data, dtype="<f8").copy()
    return arr


# =========================
# 读入数据
# =========================
E = load_dat_file(DATA_DIR / "en.dat")
I = load_dat_file(DATA_DIR / "flux.dat")
Q = load_dat_file(DATA_DIR / "qflux.dat")
U = load_dat_file(DATA_DIR / "uflux.dat")

# 如果存在就读取，用于对照
poldeg_file = DATA_DIR / "poldeg.dat"
polang_file = DATA_DIR / "polang.dat"

P_file = load_dat_file(poldeg_file) if poldeg_file.exists() else None
A_file = load_dat_file(polang_file) if polang_file.exists() else None

# 长度检查
n = len(E)
for name, arr in [("flux", I), ("qflux", Q), ("uflux", U)]:
    if len(arr) != n:
        raise ValueError(f"{name} 长度 {len(arr)} 与 en.dat 长度 {n} 不一致")

if P_file is not None and len(P_file) != n:
    raise ValueError("poldeg.dat 长度与 en.dat 不一致")
if A_file is not None and len(A_file) != n:
    raise ValueError("polang.dat 长度与 en.dat 不一致")


# =========================
# 原始点计算偏振
# =========================
P_raw = np.full_like(I, np.nan)
mask_I_pos = I > 0
P_raw[mask_I_pos] = np.sqrt(Q[mask_I_pos] ** 2 + U[mask_I_pos] ** 2) / I[mask_I_pos]
psi_raw = 0.5 * np.arctan2(U, Q)   # rad


# =========================
# 先按通量阈值选取有效数据范围
# =========================
finite_mask = np.isfinite(E) & np.isfinite(I) & np.isfinite(Q) & np.isfinite(U)
finite_mask &= (E > 0)

if not np.any(finite_mask):
    raise ValueError("没有可用的有限数据点。")

E0 = E[finite_mask]
I0 = I[finite_mask]
Q0 = Q[finite_mask]
U0 = U[finite_mask]

Imax = np.max(I0)
base_mask = I0 > (FLUX_FRAC * Imax)

E0 = E0[base_mask]
I0 = I0[base_mask]
Q0 = Q0[base_mask]
U0 = U0[base_mask]

if len(E0) == 0:
    raise ValueError("应用通量阈值后没有剩余数据，请调小 FLUX_FRAC。")


# =========================
# log 能量 bin：先 bin I,Q,U，再算 P 和 psi
# =========================
emin = E0.min()
emax = E0.max()

bins = np.logspace(np.log10(emin), np.log10(emax), NBIN + 1)

Ebin = []
Ibin = []
Qbin = []
Ubin = []

for i in range(NBIN):
    if i < NBIN - 1:
        m = (E0 >= bins[i]) & (E0 < bins[i + 1])
    else:
        m = (E0 >= bins[i]) & (E0 <= bins[i + 1])

    if np.sum(m) == 0:
        continue

    # 几何中心更适合 log 图
    e_center = np.sqrt(bins[i] * bins[i + 1])

    Ebin.append(e_center)
    Ibin.append(np.sum(I0[m]))
    Qbin.append(np.sum(Q0[m]))
    Ubin.append(np.sum(U0[m]))

Ebin = np.asarray(Ebin, dtype=float)
Ibin = np.asarray(Ibin, dtype=float)
Qbin = np.asarray(Qbin, dtype=float)
Ubin = np.asarray(Ubin, dtype=float)

if len(Ebin) == 0:
    raise ValueError("bin 之后没有数据，请减少 NBIN 或检查输入数据。")

Pbin = np.full_like(Ibin, np.nan)
pos = Ibin > 0
Pbin[pos] = np.sqrt(Qbin[pos] ** 2 + Ubin[pos] ** 2) / Ibin[pos]
psibin = 0.5 * np.arctan2(Ubin, Qbin)   # rad
psibin_deg = np.degrees(psibin)


# =========================
# 一致性检查
# =========================
print("=" * 70)
print("原始数组长度:", n)
print("有效点数:", len(E0))
print("bin 后点数:", len(Ebin))
print("E range:", emin, "to", emax)
print("I max:", np.max(I0))
print("Pbin min/max:", np.nanmin(Pbin), np.nanmax(Pbin))
print("Psi bin range (deg):", np.nanmin(psibin_deg), np.nanmax(psibin_deg))

if P_file is not None:
    diff = np.nanmax(np.abs(P_raw[mask_I_pos] - P_file[mask_I_pos]))
    print("原始 poldeg.dat 与 sqrt(Q^2+U^2)/I 最大差值 =", diff)

print("=" * 70)


# =========================
# 输出平滑后的数据
# =========================
np.savetxt(
    OUT_DIR / f"{OUT_PREFIX}_binned_stokes.txt",
    np.column_stack([Ebin, Ibin, Qbin, Ubin, Pbin, psibin_deg]),
    header="E_keV I_bin Q_bin U_bin P_bin psi_deg"
)


# =========================
# 图 1：原始 Stokes
# =========================
plt.figure(figsize=(7, 5))
plt.loglog(E0, I0, label="I")
plt.loglog(E0, np.abs(Q0), label="|Q|")
plt.loglog(E0, np.abs(U0), label="|U|")
plt.xlabel("Energy (keV)")
plt.ylabel("Flux / Stokes")
plt.title("Raw spectrum and Stokes parameters")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DIR / f"{OUT_PREFIX}_raw_stokes.png", dpi=300)
plt.close()


# =========================
# 图 2：平滑后 Stokes
# =========================
plt.figure(figsize=(7, 5))
plt.loglog(Ebin, Ibin, label="I (binned)")
plt.loglog(Ebin, np.abs(Qbin), label="|Q| (binned)")
plt.loglog(Ebin, np.abs(Ubin), label="|U| (binned)")
plt.xlabel("Energy (keV)")
plt.ylabel("Binned Flux / Stokes")
plt.title("Binned spectrum and Stokes parameters")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DIR / f"{OUT_PREFIX}_binned_stokes.png", dpi=300)
plt.close()


# =========================
# 图 3：原始偏振度 vs 平滑后偏振度
# =========================
plt.figure(figsize=(7, 5))
plt.semilogx(E0, P_raw[finite_mask][base_mask], alpha=0.35, label="raw")
plt.semilogx(Ebin, Pbin, lw=2, label="binned")
plt.xlabel("Energy (keV)")
plt.ylabel("Polarization degree")
plt.title("Polarization degree")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DIR / f"{OUT_PREFIX}_poldeg.png", dpi=300)
plt.close()


# =========================
# 图 4：偏振角
# =========================
plt.figure(figsize=(7, 5))
plt.semilogx(Ebin, psibin_deg, lw=2)
plt.xlabel("Energy (keV)")
plt.ylabel("Polarization angle (deg)")
plt.title("Binned polarization angle")
plt.tight_layout()
plt.savefig(OUT_DIR / f"{OUT_PREFIX}_polang.png", dpi=300)
plt.close()


# =========================
# 图 5：论文风格两面板
# =========================
fig = plt.figure(figsize=(7, 7))

ax1 = fig.add_subplot(2, 1, 1)
ax1.loglog(Ebin, Ebin * Ibin)
ax1.set_ylabel(r"$E \times I(E)$")
ax1.set_title("Warm corona slab polarization (binned)")

ax2 = fig.add_subplot(2, 1, 2)
ax2.semilogx(Ebin, Pbin)
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel("Polarization degree")

plt.tight_layout()
plt.savefig(OUT_DIR / f"{OUT_PREFIX}_paperstyle.png", dpi=300)
plt.close()


print("输出完成：")
print(OUT_DIR / f"{OUT_PREFIX}_binned_stokes.txt")
print(OUT_DIR / f"{OUT_PREFIX}_raw_stokes.png")
print(OUT_DIR / f"{OUT_PREFIX}_binned_stokes.png")
print(OUT_DIR / f"{OUT_PREFIX}_poldeg.png")
print(OUT_DIR / f"{OUT_PREFIX}_polang.png")
print(OUT_DIR / f"{OUT_PREFIX}_paperstyle.png")