#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一版：slab + calspec 批处理（Linux/WSL 内运行）
- 阻塞执行，无 busy-wait
- Python 修改 params.txt（不依赖 sed）
- 并行执行（ProcessPoolExecutor），workers=SLURM_CPUS_PER_TASK 或默认 40
"""

import os
import time
import random
import shutil
import subprocess
from pathlib import Path
from tqdm import tqdm

# ====== 可配置区 ======
SLAB_BIN     = "/public5/home/t6s009182/data/monk/monk_for_rhy/bin/slab"
CALSPEC_BIN  = "/public5/home/t6s009182/data/monk/monk_for_rhy/bin/calspec"

# 日志与数据根（务必用 Linux 路径）
BASE_FILE_PATH        = "/public5/home/t6s009182/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/data"
BASE_FILE_PATH_LINUX  = "/public5/home/t6s009182/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/"
PARAMS_TXT_PATH       = os.path.join(BASE_FILE_PATH_LINUX, "params.txt")

# ====== 运行外部程序（阻塞 + 抛异常） ======
def run_slab_program(slab_folder: str) -> None:
    """
    在 slab_folder 下执行 slab，失败抛 CalledProcessError
    """
    subprocess.run([SLAB_BIN], cwd=slab_folder, check=True)

def run_calspec_program(calspec_folder: str, calspec_parameter: str) -> None:
    """
    在 calspec_folder 下执行：
      calspec ../slab/ <a> <b> <c>
    参数从 calspec_parameter（形如 "A B C"）解析
    """
    args = str(calspec_parameter).split()
    cmd = [CALSPEC_BIN, "../slab/"] + args
    subprocess.run(cmd, cwd=calspec_folder, check=True)

# ====== 选点策略 ======
def choose_points(selection_method, selection_ratio, te_parts, tau_parts,
                  te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound,
                  random_seed=None):
    selected_points = []

    if selection_method == 'random':
        total_points = te_parts * tau_parts
        selected_points_count = round(total_points * selection_ratio)
        if random_seed is not None:
            random.seed(random_seed)
        all_points = [(i, j) for i in range(te_parts) for j in range(tau_parts)]
        selected_points = random.sample(all_points, selected_points_count)

    elif selection_method == 'row':
        selected_rows = round(te_parts * selection_ratio)
        selected_rows = int(selected_rows + 0.5) if (selected_rows + 0.5) > selected_rows else int(selected_rows)
        for i in range(selected_rows):
            for j in range(tau_parts):
                selected_points.append((i, j))

    elif selection_method == 'column':
        selected_columns = round(tau_parts * selection_ratio)
        selected_columns = int(selected_columns + 0.5) if (selected_columns + 0.5) > selected_columns else int(selected_columns)
        for i in range(te_parts):
            for j in range(selected_columns):
                selected_points.append((i, j))

    elif selection_method == 'all':
        for i in range(te_parts):
            for j in range(tau_parts):
                selected_points.append((i, j))

    # 边界加固
    for i in range(te_parts):
        selected_points.append((i, 0))
        selected_points.append((i, tau_parts - 1))
    for j in range(tau_parts):
        selected_points.append((0, j))
        selected_points.append((te_parts - 1, j))

    return selected_points

# ====== 文件/目录与 params 工具 ======
def _safe_make_dirs(*paths: str) -> None:
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)

def _write_params_with_updates(src_params: str, dst_params: str, te_value: float, tau_value: float) -> None:
    """
    复制 params.txt 并把第2/3行改为：
      te = <te_value>
      tau = <tau_value>
    （从1开始计行）
    """
    with open(src_params, "r", encoding="utf-8") as f:
        lines = f.readlines()
    while len(lines) < 3:
        lines.append("\n")

    lines[1] = f"te = {te_value:.3f}\n"
    lines[2] = f"tau = {tau_value:.3f}\n"

    with open(dst_params, "w", encoding="utf-8") as f:
        f.writelines(lines)

# ====== 单点任务 ======
def process_point(te_index, tau_index, te_value, tau_value, te_parts, tau_parts,
                  base_file_path_linux, params_path, calspec_parameter,
                  selected_points):

    folder_name   = f"te_{te_value:.3f}_tau_{tau_value:.3f}"
    folder        = os.path.join(base_file_path_linux, "data", folder_name)
    slab_folder   = os.path.join(folder, "slab")
    calspec_folder= os.path.join(folder, "calspec")
    dst_params    = os.path.join(slab_folder, "params.txt")

    # 目录
    _safe_make_dirs(folder, slab_folder, calspec_folder)
    # 写 params（拷贝+改第2/3行）
    _write_params_with_updates(params_path, dst_params, te_value, tau_value)

    # 跑 slab 与 calspec（阻塞，失败抛异常）
    run_slab_program(slab_folder)
    run_calspec_program(calspec_folder, calspec_parameter)

    # 若在选点集合里，返回记录路径
    if (te_index, tau_index) in selected_points:
        flux_path = os.path.join(calspec_folder, "flux.dat")
        en_path   = os.path.join(calspec_folder, "en.dat")
        de_path   = os.path.join(calspec_folder, "de.dat")
        tag       = folder_name
        return (flux_path, en_path, de_path, tag)

    return None

# ====== 主流程 ======
def create_folders(te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound,
                   te_parts, tau_parts, calspec_parameter, selection_method,
                   selection_ratio, random_seed):
    start_time = time.time()

    te_values  = [te_lower_bound  + (te_upper_bound  - te_lower_bound)  * i / (te_parts  - 1) for i in range(te_parts)]
    tau_values = [tau_lower_bound + (tau_upper_bound - tau_lower_bound) * i / (tau_parts - 1) for i in range(tau_parts)]
    num_folders = te_parts * tau_parts

    selected_points = choose_points(selection_method, selection_ratio, te_parts, tau_parts,
                                    te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound, random_seed)

    # 日志
    Path(BASE_FILE_PATH).mkdir(parents=True, exist_ok=True)
    log_file_path       = os.path.join(BASE_FILE_PATH, "control.log")
    chose_log_file_path = os.path.join(BASE_FILE_PATH, "chose.log")
    xspec_flux_log_path = os.path.join(BASE_FILE_PATH, "Xspec_flux.log")

    with open(log_file_path, "w", encoding="utf-8") as log_file, \
         open(chose_log_file_path, "w", encoding="utf-8") as chose_log_file, \
         open(xspec_flux_log_path, "w", encoding="utf-8") as xspec_flux_log_file:

        log_file.write(f"te_lower_bound: {te_lower_bound}\n")
        log_file.write(f"te_upper_bound: {te_upper_bound}\n")
        log_file.write(f"tau_lower_bound: {tau_lower_bound}\n")
        log_file.write(f"tau_upper_bound: {tau_upper_bound}\n")
        log_file.write(f"te_parts: {te_parts}\n")
        log_file.write(f"tau_parts: {tau_parts}\n")
        log_file.write(f"calspec_parameter: {calspec_parameter}\n")
        log_file.write(f"Total_number_of_folders: {num_folders}\n")

        chose_log_file.write(f"Selection Method: {selection_method}\n")
        chose_log_file.write(f"Selection Ratio: {selection_ratio}\n")
        chose_log_file.write(f"Random Seed: {random_seed}\n")

        # 并发度：Slurm 给就用，没有就默认 40
        workers = int(os.getenv("SLURM_CPUS_PER_TASK", "9"))

        from concurrent.futures import ProcessPoolExecutor, as_completed
        tasks = []
        with ProcessPoolExecutor(max_workers=workers) as executor:
            for i in range(num_folders):
                te_index  = i // tau_parts
                tau_index = i % tau_parts
                te_value  = te_values[te_index]
                tau_value = tau_values[tau_index]

                fut = executor.submit(
                    process_point, te_index, tau_index, te_value, tau_value, te_parts, tau_parts,
                    BASE_FILE_PATH_LINUX, PARAMS_TXT_PATH, calspec_parameter, selected_points
                )
                tasks.append(fut)

            for fut in tqdm(as_completed(tasks), total=len(tasks), desc="进度", unit="点"):
                try:
                    result = fut.result()
                except Exception as e:
                    # 出错写日志，不中断全局
                    log_file.write(f"[ERROR] {e}\n")
                    continue

                if result is not None:
                    flux_path, en_path, de_path, tag = result
                    xspec_flux_log_file.write(flux_path + '\n')
                    xspec_flux_log_file.write(en_path   + '\n')
                    xspec_flux_log_file.write(de_path   + '\n\n')
                    chose_log_file.write(tag + '\n')

        elapsed = time.time() - start_time
        log_file.write(f"\n运行时间: {elapsed:.3f} 秒\n")
        print(f"运行时间：{elapsed:.3f} s")

# ====== CLI ======
import argparse

def main():
    parser = argparse.ArgumentParser(description="控制 slab + calspec 批处理参数")
    parser.add_argument('--test', type=int, choices=[0, 1], required=True, help='是否测试阶段（1是，0否）')

    # 非测试参数
    parser.add_argument('--te-lower', type=float, help='te区间下限')
    parser.add_argument('--te-upper', type=float, help='te区间上限')
    parser.add_argument('--tau-lower', type=float, help='tau区间下限')
    parser.add_argument('--tau-upper', type=float, help='tau区间上限')
    parser.add_argument('--te-parts', type=int, help='te等分份数')
    parser.add_argument('--tau-parts', type=int, help='tau等分份数')
    parser.add_argument('--calspec', nargs=3, type=float, help='calspec参数（三个数字，空格分隔）')

    # 选点
    parser.add_argument('--select', choices=['row', 'column', 'random', 'all'], required=True,
                        help='选点方式：row/column/random/all')
    parser.add_argument('--ratio', type=float, default=1.0, help='数据点占比（0~1），仅 row/column/random 有效')
    parser.add_argument('--seed', type=int, default=None, help='随机种子，仅 random 有效')

    args = parser.parse_args()

    if args.test == 1:
        te_lower_bound  = 0.1
        te_upper_bound  = 0.2
        tau_lower_bound = 0.01
        tau_upper_bound = 0.02
        te_parts  = 2
        tau_parts = 2
        calspec_parameter = f"{-200} {1e-2} {1e2}"
    else:
        te_lower_bound  = args.te_lower
        te_upper_bound  = args.te_upper
        tau_lower_bound = args.tau_lower
        tau_upper_bound = args.tau_upper
        te_parts  = args.te_parts
        tau_parts = args.tau_parts
        calspec_parameter = f"{args.calspec[0]} {args.calspec[1]} {args.calspec[2]}"

    selection_method = args.select
    selection_ratio  = args.ratio if selection_method in ['row', 'column', 'random'] else 1
    random_seed      = args.seed  if selection_method == 'random' else None

    create_folders(te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound,
                   te_parts, tau_parts, calspec_parameter,
                   selection_method, selection_ratio, random_seed)

if __name__ == '__main__':
    main()


