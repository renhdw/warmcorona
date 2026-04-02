#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成两个 XSPEC atable (.mod)：
  1) warmcom_0.1-1.0_5-25_smoothed_tv_clean_Tom_new_107.mod   <- 使用 flux_clean_smoothed_tv.dat
  2) warmcom_0.1-1.0_5-25_smoothed_tv_Tom_new_107.mod         <- 使用 flux_smoothed_tv.dat (raw)
"""
import os
import re
import struct
import numpy as np
from typing import List, Tuple
from heasp import table, tableParameter, tableSpectrum

# ================= 用户配置 =================
base_dir = r"/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"
energy_file_rel = ("te_0.100_tau_10.000", "calspec", "Xspec_en.dat")  # 任一目录里的能量边界文件（长度=Nbins+1）
log_file_rel    = ("chose.log",)

# 输出文件（完整路径，精确名称）
tablefile_clean = r"/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_smoothed_tv_clean_Tom_slab_107.mod"
tablefile_raw   = r"/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_smoothed_tv_Tom_slab_107.mod"
tablefile_pure  = r"/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_pure_Tom_slab_107.mod"   # <- 新增


# — 输入光谱单位设定 —
INPUT_IS_PER_KEV = True   # 如果 flux 已经是 per MeV，设为 False

# 插值/参数设置
PARAM_INTERP_METHOD = 0   # 0=线性, 1=对数
DELTA_te  = 0.01
DELTA_tau = 0.01
# ======================================================

# 变体定义： (name, flux_filename, output_mod_path)

VARIANTS = [
    ("clean", "flux_clean_smoothed_tv.dat", tablefile_clean),
    ("raw",   "flux_smoothed_tv.dat",       tablefile_raw),
    ("pure",  "flux.dat",                   tablefile_pure),   # <- 新增：完全未加工的光谱
]


def load_dat_file(file_path: str) -> np.ndarray:
    with open(file_path, "rb") as f:
        data = f.read()
    if len(data) % 8 != 0:
        raise ValueError(f"{file_path}: 字节长度不是8的整数倍，无法按 double 解析。")
    arr = struct.unpack("<" + "d" * (len(data) // 8), data)
    return np.asarray(arr, dtype=float)

def get_parameters_from_log(log_file_path: str) -> List[Tuple[float, float]]:
    pairs: List[Tuple[float,float]] = []
    with open(log_file_path, "r") as f:
        lines = f.readlines()
    for line in lines[3:]:
        token = line.strip()
        if not token:
            continue
        m = re.match(r"^te_([0-9.]+)_tau_([0-9.]+)$", token)
        if not m:
            # 兼容带冒号的形式
            if ":" in token:
                token = token.split(":")[-1].strip()
                m = re.match(r"^te_([0-9.]+)_tau_([0-9.]+)$", token)
            if not m:
                print(f"Ignoring invalid parameter value: {line.strip()}")
                continue
        te  = float(m.group(1))
        tau = float(m.group(2))
        pairs.append((te, tau))
    return pairs

def build_table_for_variant(variant_name: str, flux_filename: str, e_edges: np.ndarray,
                            pair_list: List[Tuple[float, float]], tablefile_out: str) -> None:
    nbins = len(e_edges) - 1
    te_vals  = sorted({te for te, _ in pair_list})
    tau_vals = sorted({tau for _, tau in pair_list})

    tbl = table()
    tbl.setModelName(f"warmcom_{variant_name}")
    tbl.setModelUnits("ph/cm^2/s/MeV")   # 使用 /MeV 作为中间单位
    tbl.setEnergyUnits("keV")
    tbl.setisRedshift(True)
    tbl.setisAdditive(True)
    tbl.setisError(False)
    tbl.setEnergies(e_edges)
    tbl.setNumIntParams(2)
    tbl.setNumAddParams(0)

    # te 参数
    p_te = tableParameter()
    p_te.setName("te")
    p_te.setInterpolationMethod(PARAM_INTERP_METHOD)
    p_te.setInitialValue(float(te_vals[0]))
    p_te.setDelta(float(DELTA_te))
    p_te.setMinimum(float(min(te_vals))); p_te.setBottom(float(min(te_vals)))
    p_te.setTop(float(max(te_vals)));     p_te.setMaximum(float(max(te_vals)))
    p_te.setTabulatedValues(np.asarray(te_vals, dtype=float))
    tbl.pushParameter(p_te)

    # tau 参数
    p_tau = tableParameter()
    p_tau.setName("tau")
    p_tau.setInterpolationMethod(PARAM_INTERP_METHOD)
    p_tau.setInitialValue(float(tau_vals[0]))
    p_tau.setDelta(float(DELTA_tau))
    p_tau.setMinimum(float(min(tau_vals))); p_tau.setBottom(float(min(tau_vals)))
    p_tau.setTop(float(max(tau_vals)));     p_tau.setMaximum(float(max(tau_vals)))
    p_tau.setTabulatedValues(np.asarray(tau_vals, dtype=float))
    tbl.pushParameter(p_tau)

    miss = 0
    for te, tau in pair_list:
        folder = f"te_{te:.3f}_tau_{tau:.3f}"
        flux_path = os.path.join(base_dir, folder, "calspec", flux_filename)
        if not os.path.exists(flux_path):
            print(f"⚠️ [{variant_name}] 缺少光谱文件：{flux_path}")
            miss += 1
            continue

        flux_density = load_dat_file(flux_path)
        if len(flux_density) != nbins:
            raise ValueError(f"{flux_path}: 长度 {len(flux_density)} 与能量 bins {nbins} 不一致。")

        if INPUT_IS_PER_KEV:
            flux_density = flux_density * 1000.0  # /keV -> /MeV

        flux_density = np.asarray(flux_density, dtype=float)
        flux_density[~np.isfinite(flux_density)] = 0.0
        flux_density[flux_density < 0] = 0.0

        spec = tableSpectrum()
        spec.setParameterValues(np.array([te, tau], dtype=float))
        spec.setFlux(flux_density)
        tbl.pushSpectrum(spec)

    if miss:
        print(f"提示：变体 [{variant_name}] 中共有 {miss} 个 (te,tau) 组合缺少光谱文件。")

    ret = tbl.convertUnits()
    if ret != 0:
        raise RuntimeError(f"[{variant_name}] convertUnits() 失败，返回码 {ret}；请确认 setModelUnits 是否正确。")

    if os.path.exists(tablefile_out):
        os.remove(tablefile_out)
    status = tbl.write(tablefile_out)
    if status != 0:
        raise RuntimeError(f"[{variant_name}] 写出 {tablefile_out} 失败，status={status}")

    print(f"[{variant_name}] 已保存到: {tablefile_out}")

def main():
    energy_file_path = os.path.join(base_dir, *energy_file_rel)
    e_edges = load_dat_file(energy_file_path)
    if np.any(np.diff(e_edges) <= 0):
        raise ValueError("能量边界必须严格递增。")
    nbins = len(e_edges) - 1
    if nbins <= 0:
        raise ValueError("能量 bin 数必须 > 0。")

    log_file_path = os.path.join(base_dir, *log_file_rel)
    pair_list = get_parameters_from_log(log_file_path)
    if not pair_list:
        raise RuntimeError(f"未在日志中解析到任何 te/tau 组合：{log_file_path}")

    for variant_name, flux_filename, out_path in VARIANTS:
        print(f"开始构建变体 [{variant_name}]，flux 文件 = {flux_filename}，输出 = {out_path}")
        build_table_for_variant(variant_name, flux_filename, e_edges, pair_list, out_path)

    print("全部完成。")

if __name__ == "__main__":
    main()
