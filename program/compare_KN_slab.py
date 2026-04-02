#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import struct
import subprocess
import numpy as np

import matplotlib
matplotlib.use("Agg")   # 关键：禁止调用 Qt/X11
import matplotlib.pyplot as plt


# =========================
# 基本配置
# =========================
OUT_DIR = os.path.expanduser("~/data/monk/plot/warmcorona/test/2026.03.11/compare_KN/slab")

slab_BIN = "/home/hdw/data/monk/monk_for_rhy/bin/slab"
CALSPEC_BIN = "/home/hdw/data/monk/monk_for_rhy/bin/calspec"

BASE_PARAMS_FILE = os.path.join(OUT_DIR, "params.txt")
CALSPEC_PARAMETER = "-1000 1e-2 20"

PARAM_SETS = [
    (0.4, 20.0),
    (0.3, 16.0),
    (1.0, 10.0),
]


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


def update_params_file(base_file, out_file, te, tau, stype):
    with open(base_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    new_lines = []
    in_physical = False
    in_option = False

    found_te = False
    found_tau = False
    found_stype = False

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
        elif in_option and s.startswith("stype"):
            new_lines.append(f"stype = {stype}\n")
            found_stype = True
        else:
            new_lines.append(line)

    if not found_te:
        raise ValueError("base_params.txt 中没有找到 te = ...")
    if not found_tau:
        raise ValueError("base_params.txt 中没有找到 tau = ...")
    if not found_stype:
        raise ValueError("base_params.txt 中没有找到 stype = ...")

    with open(out_file, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


def run_one_case(case_dir, te, tau, stype):
    slab_dir = os.path.join(case_dir, "slab")
    calspec_dir = os.path.join(case_dir, "calspec")
    os.makedirs(slab_dir, exist_ok=True)
    os.makedirs(calspec_dir, exist_ok=True)

    params_out = os.path.join(slab_dir, "params.txt")
    update_params_file(BASE_PARAMS_FILE, params_out, te, tau, stype)

    run_command(slab_BIN, cwd=slab_dir)
    print(f"[OK] slab finished: te={te}, tau={tau}, stype={stype}")

    run_command(f"{CALSPEC_BIN} ../slab/ {CALSPEC_PARAMETER}", cwd=calspec_dir)
    print(f"[OK] calspec finished: te={te}, tau={tau}, stype={stype}")

    en = load_dat_file(os.path.join(calspec_dir, "en.dat"))
    flux = load_dat_file(os.path.join(calspec_dir, "flux.dat"))

    if len(en) != len(flux):
        raise ValueError(f"en.dat 和 flux.dat 长度不一致: {len(en)} vs {len(flux)}")

    return en, flux


def plot_compare(out_png, te, tau, en1, flux1, en2, flux2):
    n = min(len(en1), len(flux1), len(en2), len(flux2))
    en1 = en1[:n]
    flux1 = flux1[:n]
    en2 = en2[:n]
    flux2 = flux2[:n]

    # 检查能量网格是否一致
    if not np.allclose(en1, en2, rtol=0.0, atol=0.0):
        raise ValueError("stype=1 和 stype=2 的 en.dat 不一致，无法直接做 ratio。")

    ratio = np.full(n, np.nan)
    mask = (flux1 > 0) & np.isfinite(flux1) & np.isfinite(flux2)
    ratio[mask] = flux2[mask] / flux1[mask]

    fig = plt.figure(figsize=(8, 7))

    ax1 = plt.subplot(2, 1, 1)
    ax1.loglog(en1, flux1, label="stype=1")
    ax1.loglog(en2, flux2, label="stype=2")
    ax1.set_ylabel("Flux")
    ax1.set_title(f"te={te:.3f} keV, tau={tau:.3f}")
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)

    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    ax2.semilogx(en1, ratio)
    ax2.axhline(1.0, linestyle="--")
    ax2.set_xlabel("Energy (keV)")
    ax2.set_ylabel("stype2 / stype1")
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

        case1_dir = os.path.join(group_dir, "stype_1")
        case2_dir = os.path.join(group_dir, "stype_2")

        en1, flux1 = run_one_case(case1_dir, te, tau, 1)
        en2, flux2 = run_one_case(case2_dir, te, tau, 2)

        out_png = os.path.join(group_dir, "compare_stype1_vs_stype2.png")
        plot_compare(out_png, te, tau, en1, flux1, en2, flux2)

    print("\n全部完成。")


if __name__ == "__main__":
    main()