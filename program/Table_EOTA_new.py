

import os
import re
import json
import math
import logging
import numpy as np
import pandas as pd
from typing import List, Tuple, Optional
from heasp import table, tableParameter, tableSpectrum

# =============== 配置 ===============
BASE_DIR   = os.path.expanduser("~/data/naoc/EOTA")

# 参数表（包含 No 和各物理参数列；不含 logD/theta）
PARA_CSV  = os.path.join(BASE_DIR, "disk_para_sorted.csv")

# 存放多个 D/theta 扁平目录的根目录
SPECS_ROOT = os.path.join(BASE_DIR, "flux_grids_flat")  # 例如包含 Dlog_1.0_theta_00deg/ 等

# 输出
LOG_PATH  = os.path.join(BASE_DIR, "disk_build_table.log")
OUT_MOD   = os.path.join(BASE_DIR, "disk_EOTA.mod")
PARAM_TXT = os.path.join(BASE_DIR, "disk_table_param_report.txt")

# 能量解释：首行是能量中心（True）还是能量边界（False）
ENERGY_IS_CENTER = True

# 限制与安全
MAX_TOTAL_SPECTRA = 200000  # 输出的 tableSpectrum 条数上限
PRINT_FIRST_N     = 5       # 控制台仅打印前 N 条谱的参数与前 20 点 flux
# ====================================

DIR_REGEX = re.compile(r"^Dlog_(?P<logd>\d+(?:\.\d+)?)_theta_(?P<theta>\d{2})deg$")

def setup_logger(log_path: str) -> logging.Logger:
    logger = logging.getLogger("build_table_logger")
    logger.setLevel(logging.DEBUG)
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fh = logging.FileHandler(log_path, mode='w', encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(fmt); ch.setFormatter(fmt)
    logger.addHandler(fh); logger.addHandler(ch)
    return logger

def centers_to_edges_log(energy_centers: np.ndarray) -> np.ndarray:
    """
    将对数均匀的能量中心近似转换为边界数组（长度+1）。
    """
    logc = np.log10(energy_centers)
    widths = np.diff(logc)
    mean_w = np.mean(widths)
    loge = np.linspace(logc[0] - mean_w/2, logc[-1] + mean_w/2, len(energy_centers)+1)
    return 10**loge

def find_spec_dirs(root: str) -> List[Tuple[float, float, str]]:
    """
    扫描 root 下一级目录，匹配 Dlog_x.x_theta_YYdeg，返回 (logd, theta_deg, full_path) 列表（按 logd、theta 排序）。
    """
    items = []
    for name in os.listdir(root):
        m = DIR_REGEX.match(name)
        if m:
            logd = float(m.group("logd"))
            theta = float(m.group("theta"))  # 存为 float，heasp 参数也用 float
            items.append((logd, theta, os.path.join(root, name)))
    items.sort(key=lambda x: (x[0], x[1]))
    return items

def read_spec_csv(path: str) -> pd.DataFrame:
    """
    读取单个 spec.csv（无表头）。首行包括能量；第 0 列是 No。
    """
    return pd.read_csv(path, header=None, sep=None, engine="python")

def parse_energy_from_df(df: pd.DataFrame, energy_is_center: bool, logger: logging.Logger) -> np.ndarray:
    """
    从 DataFrame 第 1 行（索引 0）解析能量；返回边界数组（长度 = 有效能量列 + 1）。
    """
    energy_cols = df.columns[1:]  # 跳过 No
    energies = np.array([float(df.iloc[0, c]) for c in energy_cols], dtype=float)
    if not np.all(np.diff(energies) > 0):
        logger.warning("能量行不严格单调递增，后续可能报错。")
    if energy_is_center:
        edges = centers_to_edges_log(energies)
    else:
        edges = energies
    if not np.all(np.diff(edges) > 0):
        logger.warning("能量边界不严格单调递增。")
    return edges

def build_param_object_from_df_column(tbl: table, para_df: pd.DataFrame, col: str, logger: logging.Logger) -> tableParameter:
    """
    根据 para_df 的某一列（原有物理参数）构建一个 tableParameter，并 push 到 tbl。
    """
    p = tableParameter()
    p.setName(col)
    p.setInterpolationMethod(0)  # 线性
    init = float(para_df[col].iloc[0])
    p.setInitialValue(init)
    p.setDelta(0.1)
    vals = np.sort(para_df[col].astype(float).unique())
    p.setMinimum(float(vals[0])); p.setBottom(float(vals[0]))
    p.setTop(float(vals[-1])); p.setMaximum(float(vals[-1]))
    p.setTabulatedValues(vals.astype(float))
    tbl.pushParameter(p)
    if not np.all(np.diff(vals) > 0):
        logger.warning(f"参数 {col} 的唯一值数组不是严格单调递增！")
    return p

def build_discrete_param(tbl: table, name: str, values: List[float], init_idx: int = 0, delta: float = 0.1) -> tableParameter:
    """
    构建离散取值的插值参数（logD_Mpc 或 theta_deg）。
    """
    arr = np.array(sorted(set(values)), dtype=float)
    p = tableParameter()
    p.setName(name)
    p.setInterpolationMethod(0)  # 线性（heasp 会在给定 tabulated values 上插值）
    p.setInitialValue(float(arr[min(init_idx, len(arr)-1)]))
    p.setDelta(delta)
    p.setMinimum(float(arr[0])); p.setBottom(float(arr[0]))
    p.setTop(float(arr[-1])); p.setMaximum(float(arr[-1]))
    p.setTabulatedValues(arr)
    tbl.pushParameter(p)
    return p

def main():
    logger = setup_logger(LOG_PATH)
    logger.info("开始：扫描扁平目录并构建表模型")

    # 读取参数表（包含 No 和原始物理参数）
    para_df = pd.read_csv(PARA_CSV)
    if "No" not in para_df.columns:
        raise ValueError("参数 CSV 缺少 'No' 列")
    base_param_names = [c for c in para_df.columns if c != "No"]
    logger.info(f"原始参数列: {base_param_names}")

    # 扫描 spec 目录
    spec_dirs = find_spec_dirs(SPECS_ROOT)
    if not spec_dirs:
        raise RuntimeError(f"未在 {SPECS_ROOT} 下找到任何 Dlog_x.x_theta_YYdeg 目录")
    logger.info(f"发现 {len(spec_dirs)} 个 spec 目录")

    # 先读第一个 spec.csv，确定能量边界与列数
    first_logd, first_theta, first_dir = spec_dirs[0]
    first_csv = os.path.join(first_dir, "spec.csv")
    if not os.path.exists(first_csv):
        raise FileNotFoundError(f"找不到 {first_csv}")
    df0 = read_spec_csv(first_csv)
    energy_edges = parse_energy_from_df(df0, ENERGY_IS_CENTER, logger)
    n_bins = len(energy_edges) - 1
    flux_len = df0.shape[1] - 1
    if flux_len != n_bins:
        raise ValueError(f"{first_csv}: flux 列数({flux_len}) 与能量 bin({n_bins}) 不匹配")

    # heasp 表初始化
    tbl = table()
    tbl.setModelName("EOTA_norm")
    tbl.setModelUnits("ph/cm^2/s")
    tbl.setEnergyUnits("keV")
    tbl.setisRedshift(True)
    tbl.setisAdditive(True)
    tbl.setisError(False)
    tbl.setEnergies(energy_edges)

    # 推入原有物理参数
    param_objs = []
    for pname in base_param_names:
        p = build_param_object_from_df_column(tbl, para_df, pname, logger)
        param_objs.append(p)

    # 收集所有目录的 logD/theta 值列表（用于设定 tabulated values）
    logd_values = [logd for logd, _, _ in spec_dirs]
    theta_values = [theta for _, theta, _ in spec_dirs]

    # 推入新参数：logD_Mpc、theta_deg
    p_logd  = build_discrete_param(tbl, "logD_Mpc", logd_values, init_idx=0, delta=0.1)
    p_theta = build_discrete_param(tbl, "theta_deg", theta_values, init_idx=0, delta=1.0)
    param_objs.extend([p_logd, p_theta])

    # 参数名总表（写谱时按此顺序给值）
    param_names = base_param_names + ["logD_Mpc", "theta_deg"]
    logger.info(f"最终参数顺序: {param_names}")

    # 计数与打印
    total_written = 0
    printed = 0

    # 主循环：按目录遍历，每个目录中读 spec.csv，并逐行写入 tableSpectrum
    for logd, theta_deg, dpath in spec_dirs:
        spec_csv = os.path.join(dpath, "spec.csv")
        if not os.path.exists(spec_csv):
            logger.warning(f"跳过：{spec_csv} 不存在")
            continue

        df = read_spec_csv(spec_csv)

        # 检查能量行一致性（这里只校验列数；如需更严格可比较能量数组）
        if df.shape[1] - 1 != n_bins:
            raise ValueError(f"{spec_csv}: 列数与首个文件不一致")
        # 从第 1 行开始是数据行
        for ridx in range(1, len(df)):
            no = df.iloc[ridx, 0]
            try:
                no_val = float(no)
            except Exception:
                logger.warning(f"{spec_csv} 第 {ridx} 行 No 无法解析，已跳过")
                continue

            # 对应的参数行
            pr = para_df[para_df["No"] == no_val]
            if pr.empty:
                raise ValueError(f"参数表中找不到 No={no_val}（文件 {spec_csv}）")

            base_vals = np.array([float(pr.iloc[0][p]) for p in base_param_names], dtype=float)
            flux = np.array(df.iloc[ridx, 1:], dtype=float)  # 按原样使用，不做任何缩放/归一化

            # 组装参数值：原有 + 新参数
            param_values = np.concatenate([base_vals, [float(logd), float(theta_deg)]], axis=0)

            # 打印少量示例
            if printed < PRINT_FIRST_N:
                print(f"--- No={no_val:g} | logD_Mpc={logd:.1f} | theta={theta_deg:.0f} deg ---")
                for pn, pv in zip(param_names, param_values):
                    print(f"  {pn} = {pv}")
                print("  flux前20点:", flux[:20])
                print("------------------")
                printed += 1

            spec = tableSpectrum()
            spec.setParameterValues(param_values)
            spec.setFlux(flux.astype(float))
            tbl.pushSpectrum(spec)

            total_written += 1
            if total_written > MAX_TOTAL_SPECTRA:
                raise RuntimeError(f"写入的谱条数超过上限 {MAX_TOTAL_SPECTRA}，请减少目录或限制数据量。")

    # 写出 .mod
    if os.path.exists(OUT_MOD):
        os.remove(OUT_MOD)
    status = tbl.write(OUT_MOD)
    if status != 0:
        logger.error(f"写文件失败，状态码: {status}")
    else:
        logger.info(f"成功写出文件: {OUT_MOD} | 总谱条数: {total_written}")

    # 生成参数报告（摘要版）
    with open(PARAM_TXT, "w", encoding="utf-8") as f:
        f.write("表模型参数摘要\n")
        f.write("="*60 + "\n\n")
        f.write(f"原始参数列: {base_param_names}\n")
        f.write(f"logD_Mpc 取值（共 {len(set(logd_values))} 个）: {sorted(set(logd_values))[:10]}{' ...' if len(set(logd_values))>10 else ''}\n")
        f.write(f"theta_deg 取值（共 {len(set(theta_values))} 个）: {sorted(set(theta_values))[:10]}{' ...' if len(set(theta_values))>10 else ''}\n")
        f.write(f"\n总写入谱条数: {total_written}\n")
        f.write(f"输出文件: {OUT_MOD}\n")

    print(f"✅ 完成：写入 {total_written} 条 spectra -> {OUT_MOD}")
    print(f"参数摘要：{PARAM_TXT}")

if __name__ == "__main__":
    main()
