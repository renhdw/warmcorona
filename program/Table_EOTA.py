#!/usr/bin/env python3
# 只用单条归一化版本，带打印参数和flux示例到参数txt文件
import os
import numpy as np
import pandas as pd
from heasp import table, tableParameter, tableSpectrum
import logging

def setup_logger(log_path):
    """
    创建并返回一个日志记录器，同时将日志输出到文件和控制台
    """
    logger = logging.getLogger("disk_test_logger")
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_path, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def centers_to_edges_log(energy_centers):
    """
    根据能量中心点的对数坐标，计算对应的能量边界（edges）
    使bin宽度在对数空间近似均匀
    """
    log_centers = np.log10(energy_centers)
    bin_widths = np.diff(log_centers)
    mean_bin_width = np.mean(bin_widths)
    log_edges = np.linspace(log_centers[0] - mean_bin_width / 2,
                            log_centers[-1] + mean_bin_width / 2,
                            len(energy_centers) + 1)
    edges = 10 ** log_edges
    return edges

def output_param_settings_with_flux(param_list, param_names, para_df, flux_array_list, filename, print_flux_num=20):
    """
    将模型参数信息和对应归一化后的flux写入文本文件
    param_list: heasp参数对象列表
    param_names: 参数名列表
    para_df: 参数DataFrame，必须包含"No"列
    flux_array_list: flux数组列表（原始flux，未归一化）
    filename: 输出文件路径
    print_flux_num: 每条flux打印多少点
    """
    with open(filename, "w", encoding="utf-8") as f:
        f.write("模型参数设置详细信息及单条归一化Flux前若干点\n")
        f.write("="*60 + "\n\n")

        for idx, row in para_df.iterrows():
            no = row["No"]
            f.write(f"序号 No={no}\n")
            f.write("-"*60 + "\n")

            # 写参数信息
            for i, p in enumerate(param_list):
                pname = param_names[i]
                val = float(row[pname])
                f.write(f"参数 {pname}: {val}\n")

            # 单条归一化flux
            flux = flux_array_list[idx]
            flux_max = flux.max()
            flux_norm = flux / flux_max if flux_max > 0 else flux

            # 写flux前print_flux_num点
            flux_str = ", ".join([f"{v:.6e}" for v in flux_norm[:print_flux_num]])
            f.write(f"单条归一化Flux前{print_flux_num}点: [{flux_str}]\n\n")

def main():
    base_dir = os.path.expanduser("~/data/naoc/EOTA")
    log_file = os.path.join(base_dir, "disk_test.log")
    logger = setup_logger(log_file)

    para_path = os.path.join(base_dir, "disk_para_sorted.csv")
    spec_path = os.path.join(base_dir, "disk_spec_sorted.csv")

    logger.info("开始读取参数文件和光谱文件")
    spec_df = pd.read_csv(spec_path)
    para_df = pd.read_csv(para_path)

    max_rows = 8000
    if len(spec_df) > max_rows:
        logger.info(f"原始数据行数: {len(spec_df)}，截取前{max_rows}行")
        spec_df = spec_df.iloc[:max_rows].reset_index(drop=True)
        para_df = para_df.iloc[:max_rows].reset_index(drop=True)

    energy_cols = spec_df.columns[1:]
    energy_is_center = True

    energy_vals = np.array([float(c) for c in energy_cols])
    logger.info(f"能量中心点数组大小: {energy_vals.shape}")

    if not np.all(np.diff(energy_vals) > 0):
        logger.warning("能量中心点不是严格单调递增！")

    if energy_is_center:
        energy_edges = centers_to_edges_log(energy_vals)
        logger.info("根据能量中心点计算了能量边界")
    else:
        energy_edges = energy_vals
        logger.info("直接使用能量边界数组")

    logger.info(f"能量边界数组大小: {energy_edges.shape}")
    if not np.all(np.diff(energy_edges) > 0):
        logger.warning("能量边界不是严格单调递增！")

    n_flux_bins = len(energy_edges) - 1
    flux_len = len(spec_df.iloc[0]) - 1
    logger.info(f"能量bin数: {n_flux_bins}, 光谱flux长度: {flux_len}")

    if flux_len != n_flux_bins:
        logger.error(f"flux长度 {flux_len} 与能量bin数 {n_flux_bins} 不匹配，程序退出")
        raise ValueError(f"flux长度 {flux_len} 与能量bin数 {n_flux_bins} 不匹配")

    # 预先读取所有flux到列表
    flux_array_list = []
    for idx, row in spec_df.iterrows():
        flux = np.array(row[1:], dtype=float)
        flux_array_list.append(flux)

    tbl = table()
    tbl.setModelName("EOTA_norm")
    tbl.setModelUnits("ph/cm^2/s")
    tbl.setEnergyUnits("keV")
    tbl.setisRedshift(True)
    tbl.setisAdditive(True)
    tbl.setisError(False)
    tbl.setEnergies(energy_edges)

    param_names = [c for c in para_df.columns if c != "No"]
    tbl.setNumIntParams(len(param_names))
    tbl.setNumAddParams(0)
    logger.info(f"模型参数名: {param_names}")

    param_objs = []
    for pname in param_names:
        p = tableParameter()
        p.setName(pname)
        p.setInterpolationMethod(0)
        p.setInitialValue(float(para_df[pname].iloc[0]))
        p.setDelta(0.1)
        p_min = float(para_df[pname].min())
        p_max = float(para_df[pname].max())
        p.setMinimum(p_min)
        p.setBottom(p_min)
        p.setTop(p_max)
        p.setMaximum(p_max)
        unique_vals = np.sort(para_df[pname].unique())
        p.setTabulatedValues(unique_vals)
        tbl.pushParameter(p)
        param_objs.append(p)
        logger.info(f"参数 {pname} 初始化: 初始值={p.getInitialValue()}, 范围=({p_min}, {p_max}), 唯一值数量={len(unique_vals)})")

        if not np.all(np.diff(unique_vals) > 0):
            logger.warning(f"参数 {pname} 的唯一值数组不是严格单调递增！")

    param_setting_file = os.path.join(base_dir, "disk_test_param_settings.txt")
    output_param_settings_with_flux(param_objs, param_names, para_df, flux_array_list, param_setting_file, print_flux_num=20)
    logger.info(f"参数和归一化Flux信息已写入: {param_setting_file}")

    flux_threshold = 1e-100
    print_count = 5
    printed = 0

    for idx, row in spec_df.iterrows():
        no = row["No"]
        flux = flux_array_list[idx]

        # 单条归一化
        flux_max = flux.max()
        flux_norm = flux / flux_max if flux_max > 0 else flux

        if flux_norm.max() < flux_threshold:
            logger.info(f"No={no} 归一化后最大flux={flux_norm.max()} < {flux_threshold}，数据跳过")
            continue

        para_row = para_df[para_df["No"] == no]
        if para_row.empty:
            logger.error(f"参数文件中找不到 No={no}")
            raise ValueError(f"参数文件中找不到 No={no}")
        para_vals = para_row.iloc[0]
        param_values = np.array([float(para_vals[p]) for p in param_names])
        logger.info(f"No={no} 参数值={param_values} 归一化后flux长度={len(flux_norm)}")

        for i, pname in enumerate(param_names):
            vals = tbl.getParameter(i).getTabulatedValues()
            if param_values[i] < vals[0] or param_values[i] > vals[-1]:
                logger.warning(f"No={no} 参数 {pname}={param_values[i]} 超出唯一值范围 ({vals[0]}, {vals[-1]})")

        if printed < print_count:
            print(f"--- No={no} ---")
            print("参数:")
            for pname, val in zip(param_names, param_values):
                print(f"  {pname} = {val}")
            print(f"归一化flux前20点: {flux_norm[:20]}")
            print("------------------")
            printed += 1

        spec = tableSpectrum()
        spec.setParameterValues(param_values)
        spec.setFlux(flux_norm)
        tbl.pushSpectrum(spec)

    outfile = os.path.join(base_dir, "disk_EOTA.mod")
    if os.path.exists(outfile):
        os.remove(outfile)

    status = tbl.write(outfile)
    if status != 0:
        logger.error(f"写文件失败，状态码: {status}")
    else:
        logger.info(f"成功写出文件: {outfile}")

if __name__ == "__main__":
    main()
