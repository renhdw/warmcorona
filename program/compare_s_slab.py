#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import struct
import subprocess
import numpy as np

import matplotlib
matplotlib.use("Agg")   # 禁止调用 Qt/X11
import matplotlib.pyplot as plt


# =========================
# 基本配置
# =========================
OUT_DIR = os.path.expanduser("~/data/monk/plot/warmcorona/test/2026.03.11/compare_s/slab")

slab_BIN = "/home/hdw/data/monk/monk_for_rhy/bin/slab"
CALSPEC_BIN = "/home/hdw/data/monk/monk_for_rhy/bin/calspec"

BASE_PARAMS_FILE = os.path.join(OUT_DIR, "params.txt")
CALSPEC_PARAMETER = "-1000 1e-2 20"

PARAM_SETS = [
    (0.4, 20.0),
    (0.3, 16.0),
    (1.0, 10.0),
]

# 固定 stype，只比较 s
STYPE_FIXED = 1

# 要比较的 s
S_VALUES = [1, 5, 10, 100]


# =========================
# 工具函数
# =========================
def run_command(command, cwd=None):
    result = subprocess.run(
        ["bash", "-c", command],
        cwd=cwd,
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("命令失败：")
        print(command)
        print(result.stderr)
        raise RuntimeError("command failed")
    return result


def load_dat_file(file_path: str) -> np.ndarray:
    with open(file_path, "rb") as f:
        data = f.read()
    if len(data) % 8 != 0:
        raise ValueError(f"{file_path}: 字节长度不是8的整数倍，无法按 double 解析。")
    arr = struct.unpack("<" + "d" * (len(data) // 8), data)
    return np.asarray(arr, dtype=float)


def update_params_file(base_file, out_file, te, tau, stype, s_value):
    with open(base_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    new_lines = []
    in_physical = False
    in_option = False

    found_te = False
    found_tau = False
    found_stype = False
    found_s = False

    for line in lines:
        s = line.strip()

        if s.startswith("[") and s.endswith("]"):
            sec = s.lower()
            in_physical = (sec == "[physical]")
            in_option = (sec == "[option]")
            new_lines.append(line)
            continue

        if in_physical and s.startswith("te"):
            new_lines.append(f"te = {te}\n")
            found_te = True
        elif in_physical and s.startswith("tau"):
            new_lines.append(f"tau = {tau}\n")
            found_tau = True
        elif in_physical and s.startswith("s"):
            new_lines.append(f"s = {s_value}\n")
            found_s = True
        elif in_option and s.startswith("stype"):
            new_lines.append(f"stype = {stype}\n")
            found_stype = True
        else:
            new_lines.append(line)

    if not found_te:
        raise ValueError("base params.txt 中没有找到 te = ...")
    if not found_tau:
        raise ValueError("base params.txt 中没有找到 tau = ...")
    if not found_stype:
        raise ValueError("base params.txt 中没有找到 stype = ...")
    if not found_s:
        raise ValueError("base params.txt 中没有找到 s = ...")

    with open(out_file, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


def run_one_case(case_dir, te, tau, stype, s_value):
    slab_dir = os.path.join(case_dir, "slab")
    calspec_dir = os.path.join(case_dir, "calspec")
    os.makedirs(slab_dir, exist_ok=True)
    os.makedirs(calspec_dir, exist_ok=True)

    params_out = os.path.join(slab_dir, "params.txt")
    update_params_file(BASE_PARAMS_FILE, params_out, te, tau, stype, s_value)

    run_command(slab_BIN, cwd=slab_dir)
    print(f"[OK] slab finished: te={te}, tau={tau}, stype={stype}, s={s_value}")

    run_command(f"{CALSPEC_BIN} ../slab/ {CALSPEC_PARAMETER}", cwd=calspec_dir)
    print(f"[OK] calspec finished: te={te}, tau={tau}, stype={stype}, s={s_value}")

    en = load_dat_file(os.path.join(calspec_dir, "en.dat"))
    flux = load_dat_file(os.path.join(calspec_dir, "flux.dat"))

    if len(en) != len(flux):
        raise ValueError(f"en.dat 和 flux.dat 长度不一致: {len(en)} vs {len(flux)}")

    return en, flux


def plot_compare(out_png, te, tau, results_dict):
    """
    results_dict:
        {
            1: (en, flux),
            5: (en, flux),
            10: (en, flux),
            100: (en, flux),
        }
    """
    s_ref = S_VALUES[0]
    en_ref, flux_ref = results_dict[s_ref]

    # 检查所有能量网格一致
    for s_value in S_VALUES[1:]:
        en_i, _ = results_dict[s_value]
        if not np.allclose(en_ref, en_i, rtol=0.0, atol=0.0):
            raise ValueError(f"s={s_value} 和 s={s_ref} 的 en.dat 不一致，无法直接做 ratio。")

    fig = plt.figure(figsize=(8, 8))

    ax1 = plt.subplot(2, 1, 1)
    for s_value in S_VALUES:
        en_i, flux_i = results_dict[s_value]
        ax1.loglog(en_i, flux_i, label=f"s={s_value}")
    ax1.set_ylabel("Flux")
    ax1.set_title(f"te={te:.3f} keV, tau={tau:.3f}, stype={STYPE_FIXED}")
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)

    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    for s_value in S_VALUES:
        en_i, flux_i = results_dict[s_value]
        ratio = np.full(len(en_ref), np.nan)
        mask = (flux_ref > 0) & np.isfinite(flux_ref) & np.isfinite(flux_i)
        ratio[mask] = flux_i[mask] / flux_ref[mask]
        ax2.semilogx(en_ref, ratio, label=f"s={s_value} / s={s_ref}")
    ax2.axhline(1.0, linestyle="--")
    ax2.set_xlabel("Energy (keV)")
    ax2.set_ylabel("Ratio")
    ax2.legend()
    ax2.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[OK] saved plot: {out_png}")


# =========================
# 主程序
# =========================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(BASE_PARAMS_FILE):
        raise FileNotFoundError(
            f"找不到基础参数文件：{BASE_PARAMS_FILE}\n"
            f"请先放一个 params.txt"
        )

    for te, tau in PARAM_SETS:
        group_name = f"te_{te:.3f}_tau_{tau:.3f}"
        group_dir = os.path.join(OUT_DIR, group_name)
        os.makedirs(group_dir, exist_ok=True)

        results_dict = {}

        for s_value in S_VALUES:
            case_dir = os.path.join(group_dir, f"s_{s_value}")
            en, flux = run_one_case(case_dir, te, tau, STYPE_FIXED, s_value)
            results_dict[s_value] = (en, flux)

        out_png = os.path.join(group_dir, "compare_s_values.png")
        plot_compare(out_png, te, tau, results_dict)

    print("\n全部完成。")


if __name__ == "__main__":
    main()