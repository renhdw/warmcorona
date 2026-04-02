#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
计算 nthComp 的 tau
公式：
tau = sqrt( 1 / (theta_e * ((Gamma + 0.5)^2 - 9/4)) ) - 1/3
其中 theta_e = kTe_nth / 511 keV
"""

# ==============================
# 输入参数
# ==============================
gamma    = float(input("请输入 nthComp 的 Gamma 值: "))
kTe_nth  = float(input("请输入 nthComp 的 kT_e 值 (keV): "))

# ==============================
# 计算 tau
# ==============================
me_keV = 511.0                      # 电子静能量 keV
theta_e = kTe_nth / me_keV
term = (gamma + 0.5)**2 - 2.25       # 2.25 = 9/4

if term <= 0 or theta_e <= 0:
    raise RuntimeError("参数不物理：需满足 (Gamma+0.5)^2 > 9/4 且 kT_e > 0。")

tau_nth = (1.0 / (theta_e * term))**0.5 - 1.0/3.0

print(f"[nthComp] 由 Γ={gamma:.3f}, kT_e={kTe_nth:.3f} keV 推得 tau ≈ {tau_nth:.3f}")
