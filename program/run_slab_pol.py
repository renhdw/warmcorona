#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_slab_pol.py

功能
------
给定一组 te 和 tau，自动完成以下流程：

1. 在指定数据目录下建立参数点文件夹：
      te_xxx_tau_xxx

2. 若该参数点尚未完成 slab 计算，则：
   - 从模板 params.txt 复制参数文件
   - 自动修改 te 和 tau
   - 调用 MONK 的 slab 程序

3. slab 运行后，若存在：
      qarr.dat, uarr.dat
   则自动复制为：
      qweight.dat, uweight.dat
   因为 calspec 只有在检测到 qweight.dat / uweight.dat 时，
   才会计算偏振相关输出（qflux/uflux/poldeg/polang）。

4. 运行一次“总谱版” calspec：
   输出目录：
      te_xxx_tau_xxx/calspec

5. 再按 3 度一个角度 bin，运行所有角度版 calspec：
      00-03, 03-06, ..., 87-90
   输出目录：
      te_xxx_tau_xxx/calspec_00-03
      te_xxx_tau_xxx/calspec_03-06
      ...
      te_xxx_tau_xxx/calspec_87-90

6. 【新增功能】可选地按 mu = cos(theta) 均匀分箱，运行额外一套 calspec：
      calspec_mu_00
      calspec_mu_01
      ...
      calspec_mu_{nmu-1}
   使用参数：
      --run-mu
      --nmu
      --force-mu

7. 最后读取原始角度结果，画出：
   - 不同角度的 flux(E)
   - 不同角度的 polarization degree
   - 不同角度的 polarization angle
   - 组合图（E*F(E) + P(E)）

输出图统一保存到：
   ~/data/monk/plot/warmcorona/test/2026.03.30/pol/te_xxx_tau_xxx/

使用示例
--------
python3 run_slab_pol.py --te 0.300 --tau 15.000

新增 mu 分箱：
python3 run_slab_pol.py --te 0.300 --tau 15.000 --run-mu --nmu 30

可选参数：
--force-calspec
    即使总谱 calspec 已存在，也强制重跑

--force-angle
    即使角度 calspec 已存在，也强制重跑

--force-slab
    即使 slab 结果已存在，也强制重跑 slab（慎用）

--run-mu
    额外生成 mu=cos(theta) 均匀分箱的 calspec 输出

--nmu
    mu 分箱个数，默认 30

--force-mu
    即使 mu 分箱 calspec 已存在，也强制重跑

说明
----
1. 本脚本默认：
   - slab 可执行程序路径：
       /home/hdw/data/monk/monk_for_rhy/bin/slab
   - calspec 可执行程序路径：
       /home/hdw/data/monk/monk_for_rhy/bin/calspec
   - 模板参数文件：
       /home/hdw/data/monk/plot/warmcorona/test/polarization/data/10_7/params.txt

2. 本脚本默认 calspec 参数为：
       ne   = -10000
       emin = 0.01
       emax = 600.0
   即对数能量 bin，0.01--600 keV

3. 原始角度 bin 采用：
       (0, 3], (3, 6], ..., (87, 90]
   与 calspec.cpp 中的 imin/imax 逻辑一致。

4. mu 分箱采用：
       mu in [0,1] 均匀分成 nmu 个 bin
   然后映射成 theta 区间：
       theta_min = arccos(mu_max)
       theta_max = arccos(mu_min)

5. 读取 .dat 文件时，默认按 little-endian double (<f8) 读取。
"""

import argparse
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# 基础路径设置
# ============================================================
MONK_BIN_DIR = Path("/home/hdw/data/monk/monk_for_rhy/bin")
SLAB_EXE = MONK_BIN_DIR / "slab"
CALSPEC_EXE = MONK_BIN_DIR / "calspec"

BASE_DIR = Path("/home/hdw/data/monk/plot/warmcorona/test/polarization/data/10_7")
PARAM_TEMPLATE = BASE_DIR / "params.txt"

PLOT_ROOT = Path("~/data/monk/plot/warmcorona/test/2026.04.01/pol").expanduser()


# ============================================================
# 默认 calspec 参数
# ============================================================
DEFAULT_NE = -10000
DEFAULT_EMIN = 0.01
DEFAULT_EMAX =600.0


# ============================================================
# 通用工具函数
# ============================================================
def run_bash_command(command: str, cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess:
    """
    运行 shell 命令，并可选择在失败时报错。
    """
    result = subprocess.run(
        ["bash", "-lc", command],
        cwd=str(cwd) if cwd is not None else None,
        capture_output=True,
        text=True
    )

    if check and result.returncode != 0:
        print("命令执行失败：")
        print(command)
        print("stdout:")
        print(result.stdout)
        print("stderr:")
        print(result.stderr)
        raise RuntimeError("shell 命令执行失败")

    return result


def ensure_dir(path: Path) -> None:
    """
    若目录不存在则创建。
    """
    path.mkdir(parents=True, exist_ok=True)


def copy_file(src: Path, dst: Path) -> None:
    """
    覆盖复制文件。
    """
    shutil.copy2(src, dst)


def load_dat_file(file_path: Path) -> np.ndarray:
    """
    按 little-endian double 读取 MONK/calspec 输出的二进制 .dat 文件。
    """
    raw = file_path.read_bytes()
    if len(raw) % 8 != 0:
        raise ValueError(f"{file_path} 的字节长度不是 8 的整数倍，无法按 double 解析。")
    return np.frombuffer(raw, dtype="<f8").copy()


def save_text_log(path: Path, text: str) -> None:
    """
    保存文本日志。
    """
    path.write_text(text, encoding="utf-8")


def format_tag(te: float, tau: float) -> str:
    """
    统一参数点命名格式。
    """
    return f"te_{te:.3f}_tau_{tau:.3f}"


# ============================================================
# params.txt 处理
# ============================================================
def prepare_params_file(target_dir: Path, te: float, tau: float) -> Path:
    """
    将模板 params.txt 复制到目标目录，并修改 te 与 tau。
    """
    if not PARAM_TEMPLATE.exists():
        raise FileNotFoundError(f"未找到模板参数文件：{PARAM_TEMPLATE}")

    ensure_dir(target_dir)

    out_params = target_dir / "params.txt"
    text = PARAM_TEMPLATE.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    new_lines = []
    te_done = False
    tau_done = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("te"):
            new_lines.append(f"te = {te:.3f}")
            te_done = True
        elif stripped.startswith("tau"):
            new_lines.append(f"tau = {tau:.3f}")
            tau_done = True
        else:
            new_lines.append(line)

    if not te_done:
        new_lines.append(f"te = {te:.3f}")
    if not tau_done:
        new_lines.append(f"tau = {tau:.3f}")

    out_params.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return out_params


# ============================================================
# slab 运行与检查
# ============================================================
def slab_outputs_exist(point_dir: Path) -> bool:
    """
    判断 slab 是否已经产生了关键输出。
    这里以 slab/en0.dat 和 slab/weight.dat 为最低判据。
    """
    slab_dir = point_dir / "slab"
    return (slab_dir / "en0.dat").exists() and (slab_dir / "weight.dat").exists()


def run_slab(point_dir: Path, force: bool = False) -> Path:
    """
    在 point_dir/slab 中运行 slab。
    """
    slab_dir = point_dir / "slab"
    ensure_dir(slab_dir)

    if slab_outputs_exist(point_dir) and not force:
        print(f"[跳过 slab] 已检测到 slab 结果：{slab_dir}")
        return slab_dir

    if not SLAB_EXE.exists():
        raise FileNotFoundError(f"未找到 slab 可执行文件：{SLAB_EXE}")

    print(f"[运行 slab] {slab_dir}")
    cmd = f"cd '{slab_dir}' && '{SLAB_EXE}'"
    result = run_bash_command(cmd, check=True)

    log_text = []
    log_text.append("COMMAND:")
    log_text.append(cmd)
    log_text.append("\nSTDOUT:")
    log_text.append(result.stdout)
    log_text.append("\nSTDERR:")
    log_text.append(result.stderr)
    save_text_log(slab_dir / "slab_run.log", "\n".join(log_text))

    return slab_dir


def ensure_qu_as_qweight(slab_dir: Path) -> None:
    """
    若 slab 输出的是 qarr.dat/uarr.dat，则复制成 qweight.dat/uweight.dat。
    """
    qarr = slab_dir / "qarr.dat"
    uarr = slab_dir / "uarr.dat"
    qweight = slab_dir / "qweight.dat"
    uweight = slab_dir / "uweight.dat"

    if qarr.exists():
        copy_file(qarr, qweight)
        print(f"[复制] {qarr.name} -> {qweight.name}")

    if uarr.exists():
        copy_file(uarr, uweight)
        print(f"[复制] {uarr.name} -> {uweight.name}")


# ============================================================
# calspec 运行
# ============================================================
def calspec_outputs_exist(out_dir: Path) -> bool:
    """
    判断一次 calspec 是否已经有基本输出。
    """
    return (out_dir / "en.dat").exists() and (out_dir / "flux.dat").exists() and (out_dir / "de.dat").exists()


def run_calspec_total(slab_dir: Path, point_dir: Path, ne: int, emin: float, emax: float,
                      force: bool = False) -> Path:
    """
    运行总谱版 calspec：
        calspec slab_dir ne emin emax
    输出放在：
        point_dir/calspec
    """
    if not CALSPEC_EXE.exists():
        raise FileNotFoundError(f"未找到 calspec 可执行文件：{CALSPEC_EXE}")

    out_dir = point_dir / "calspec"
    ensure_dir(out_dir)

    if calspec_outputs_exist(out_dir) and not force:
        print(f"[跳过总谱 calspec] 已存在：{out_dir}")
        return out_dir

    cmd = f"cd '{out_dir}' && '{CALSPEC_EXE}' '{slab_dir}' {ne} {emin} {emax}"
    print(f"[运行总谱 calspec] {out_dir}")
    result = run_bash_command(cmd, check=True)

    log_text = []
    log_text.append("COMMAND:")
    log_text.append(cmd)
    log_text.append("\nSTDOUT:")
    log_text.append(result.stdout)
    log_text.append("\nSTDERR:")
    log_text.append(result.stderr)
    save_text_log(out_dir / "calspec_run.log", "\n".join(log_text))

    return out_dir


def angle_dir_name(imin: int, imax: int) -> str:
    """
    原始 theta 角度目录命名，如 calspec_00-03
    """
    return f"calspec_{imin:02d}-{imax:02d}"


def run_calspec_angle(slab_dir: Path, point_dir: Path, imin: int, imax: int,
                      ne: int, emin: float, emax: float, force: bool = False) -> Path:
    """
    运行角度版 calspec：
        calspec slab_dir ne emin emax imin imax
    输出放在：
        point_dir/calspec_xx-yy
    """
    out_dir = point_dir / angle_dir_name(imin, imax)
    ensure_dir(out_dir)

    if calspec_outputs_exist(out_dir) and not force:
        print(f"[跳过角度 calspec] 已存在：{out_dir}")
        return out_dir

    cmd = f"cd '{out_dir}' && '{CALSPEC_EXE}' '{slab_dir}' {ne} {emin} {emax} {imin} {imax}"
    print(f"[运行角度 calspec] {imin:02d}-{imax:02d}")
    result = run_bash_command(cmd, check=True)

    log_text = []
    log_text.append("COMMAND:")
    log_text.append(cmd)
    log_text.append("\nSTDOUT:")
    log_text.append(result.stdout)
    log_text.append("\nSTDERR:")
    log_text.append(result.stderr)
    save_text_log(out_dir / "calspec_run.log", "\n".join(log_text))

    return out_dir


# ============================================================
# 新增：mu = cos(theta) 均匀分箱 calspec
# ============================================================
def mu_dir_name(idx: int) -> str:
    """
    mu 分箱目录命名，如 calspec_mu_00
    """
    return f"calspec_mu_{idx:02d}"


def build_mu_uniform_bins(nmu: int):
    """
    在 mu = cos(theta) ∈ [0,1] 上均匀分 nmu 个 bin。

    返回列表，每个元素是字典：
    {
        "idx": j,
        "mu_min": ...,
        "mu_max": ...,
        "theta_min": ...,
        "theta_max": ...
    }

    其中：
        theta_min = arccos(mu_max)
        theta_max = arccos(mu_min)

    因为 theta 从 face-on 的 0° 到 edge-on 的 90°，
    而 mu = cos(theta) 是递减的。
    """
    if nmu <= 0:
        raise ValueError("nmu 必须为正整数。")

    mu_edges = np.linspace(0.0, 1.0, nmu + 1)
    bins = []

    for j in range(nmu):
        mu_min = float(mu_edges[j])
        mu_max = float(mu_edges[j + 1])

        theta_min = math.degrees(math.acos(mu_max))
        theta_max = math.degrees(math.acos(mu_min))

        bins.append({
            "idx": j,
            "mu_min": mu_min,
            "mu_max": mu_max,
            "theta_min": theta_min,
            "theta_max": theta_max,
        })

    return bins


def write_mu_bin_metadata(out_dir: Path, mu_info: dict) -> None:
    """
    在 mu 分箱目录里写一个说明文件，方便之后检查。
    """
    text = (
        f"mu_index   = {mu_info['idx']}\n"
        f"mu_min     = {mu_info['mu_min']:.8f}\n"
        f"mu_max     = {mu_info['mu_max']:.8f}\n"
        f"theta_min  = {mu_info['theta_min']:.8f} deg\n"
        f"theta_max  = {mu_info['theta_max']:.8f} deg\n"
    )
    save_text_log(out_dir / "mu_bin_info.txt", text)


def run_calspec_mu_bin(slab_dir: Path, point_dir: Path, mu_info: dict,
                       ne: int, emin: float, emax: float, force: bool = False) -> Path:
    """
    运行单个 mu 均匀分箱的 calspec。

    这里 calspec 仍然只接受 theta 区间，所以将：
        mu_min, mu_max
    映射为：
        theta_min = arccos(mu_max)
        theta_max = arccos(mu_min)

    输出目录：
        point_dir/calspec_mu_xx
    """
    if not CALSPEC_EXE.exists():
        raise FileNotFoundError(f"未找到 calspec 可执行文件：{CALSPEC_EXE}")

    out_dir = point_dir / mu_dir_name(mu_info["idx"])
    ensure_dir(out_dir)

    if calspec_outputs_exist(out_dir) and not force:
        print(f"[跳过 mu 分箱 calspec] 已存在：{out_dir}")
        write_mu_bin_metadata(out_dir, mu_info)
        return out_dir

    theta_min = mu_info["theta_min"]
    theta_max = mu_info["theta_max"]

    cmd = f"cd '{out_dir}' && '{CALSPEC_EXE}' '{slab_dir}' {ne} {emin} {emax} {theta_min:.8f} {theta_max:.8f}"
    print(
        f"[运行 mu 分箱 calspec] "
        f"idx={mu_info['idx']:02d} "
        f"mu=[{mu_info['mu_min']:.5f}, {mu_info['mu_max']:.5f}] "
        f"theta=[{theta_min:.3f}, {theta_max:.3f}]"
    )
    result = run_bash_command(cmd, check=True)

    log_text = []
    log_text.append("COMMAND:")
    log_text.append(cmd)
    log_text.append("\nSTDOUT:")
    log_text.append(result.stdout)
    log_text.append("\nSTDERR:")
    log_text.append(result.stderr)
    save_text_log(out_dir / "calspec_run.log", "\n".join(log_text))

    write_mu_bin_metadata(out_dir, mu_info)
    return out_dir


def run_all_calspec_mu_bins(slab_dir: Path, point_dir: Path,
                            ne: int, emin: float, emax: float,
                            nmu: int = 30, force: bool = False) -> None:
    """
    批量运行所有 mu 均匀分箱的 calspec。
    """
    mu_bins = build_mu_uniform_bins(nmu)

    summary_lines = []
    summary_lines.append(f"nmu = {nmu}")
    summary_lines.append("idx  mu_min         mu_max         theta_min(deg)  theta_max(deg)")

    for mu_info in mu_bins:
        run_calspec_mu_bin(
            slab_dir=slab_dir,
            point_dir=point_dir,
            mu_info=mu_info,
            ne=ne,
            emin=emin,
            emax=emax,
            force=force,
        )

        summary_lines.append(
            f"{mu_info['idx']:02d}   "
            f"{mu_info['mu_min']:.8f}   "
            f"{mu_info['mu_max']:.8f}   "
            f"{mu_info['theta_min']:.8f}   "
            f"{mu_info['theta_max']:.8f}"
        )

    save_text_log(point_dir / "mu_bins_summary.txt", "\n".join(summary_lines))


# ============================================================
# 偏振计算辅助
# ============================================================
def calc_pol_from_stokes(I: np.ndarray, Q: np.ndarray, U: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    根据 I, Q, U 计算偏振度与偏振角。
    偏振角输出为 degree。
    """
    P = np.full_like(I, np.nan, dtype=float)
    psi_deg = np.full_like(I, np.nan, dtype=float)

    pos = I > 0
    polflux = np.sqrt(Q[pos] ** 2 + U[pos] ** 2)
    P[pos] = polflux / I[pos]
    psi_deg[pos] = np.degrees(0.5 * np.arctan2(U[pos], Q[pos]))

    return P, psi_deg


def read_calspec_dir(calspec_dir: Path) -> dict:
    """
    读取一个 calspec 输出目录，返回字典。
    """
    data = {
        "dir": calspec_dir,
        "en": load_dat_file(calspec_dir / "en.dat"),
        "de": load_dat_file(calspec_dir / "de.dat"),
        "flux": load_dat_file(calspec_dir / "flux.dat"),
    }

    qfile = calspec_dir / "qflux.dat"
    ufile = calspec_dir / "uflux.dat"
    pfile = calspec_dir / "poldeg.dat"
    afile = calspec_dir / "polang.dat"

    if qfile.exists() and ufile.exists():
        q = load_dat_file(qfile)
        u = load_dat_file(ufile)
        data["qflux"] = q
        data["uflux"] = u

        if pfile.exists():
            data["poldeg"] = load_dat_file(pfile)
        else:
            data["poldeg"], _ = calc_pol_from_stokes(data["flux"], q, u)

        if afile.exists():
            data["polang"] = np.degrees(load_dat_file(afile))
        else:
            _, data["polang"] = calc_pol_from_stokes(data["flux"], q, u)

    return data


# ============================================================
# 原始绘图功能（保持不变）
# ============================================================
def plot_angle_flux(angle_data: list[tuple[str, dict]], out_dir: Path, tag: str) -> None:
    """
    画不同角度的 flux(E)。
    """
    plt.figure(figsize=(7.5, 5.5))
    for label, data in angle_data:
        e = data["en"]
        f = data["flux"]
        mask = np.isfinite(e) & np.isfinite(f) & (e > 0) & (f > 0)
        if np.any(mask):
            plt.loglog(e[mask], f[mask], label=label)

    plt.xlabel("Energy (keV)")
    plt.ylabel("Flux")
    plt.title(f"Angle-resolved flux: {tag}")
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(out_dir / f"{tag}_angle_flux.png", dpi=300)
    plt.close()


def plot_angle_poldeg(angle_data: list[tuple[str, dict]], out_dir: Path, tag: str) -> None:
    """
    画不同角度的 polarization degree。
    """
    plt.figure(figsize=(7.5, 5.5))
    for label, data in angle_data:
        if "poldeg" not in data:
            continue
        e = data["en"]
        p = data["poldeg"]
        mask = np.isfinite(e) & np.isfinite(p) & (e > 0)
        if np.any(mask):
            plt.semilogx(e[mask], p[mask], label=label)

    plt.xlabel("Energy (keV)")
    plt.ylabel("Polarization degree")
    plt.title(f"Angle-resolved polarization degree: {tag}")
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(out_dir / f"{tag}_angle_poldeg.png", dpi=300)
    plt.close()


def plot_angle_polang(angle_data: list[tuple[str, dict]], out_dir: Path, tag: str) -> None:
    """
    画不同角度的 polarization angle。
    """
    plt.figure(figsize=(7.5, 5.5))
    for label, data in angle_data:
        if "polang" not in data:
            continue
        e = data["en"]
        a = data["polang"]
        mask = np.isfinite(e) & np.isfinite(a) & (e > 0)
        if np.any(mask):
            plt.semilogx(e[mask], a[mask], label=label)

    plt.xlabel("Energy (keV)")
    plt.ylabel("Polarization angle (deg)")
    plt.title(f"Angle-resolved polarization angle: {tag}")
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(out_dir / f"{tag}_angle_polang.png", dpi=300)
    plt.close()


def plot_angle_paperstyle(angle_data: list[tuple[str, dict]], out_dir: Path, tag: str) -> None:
    """
    画论文风格两面板图：
    上：E * flux(E)
    下：polarization degree
    """
    fig = plt.figure(figsize=(7.5, 7.0))

    ax1 = fig.add_subplot(2, 1, 1)
    for label, data in angle_data:
        e = data["en"]
        f = data["flux"]
        mask = np.isfinite(e) & np.isfinite(f) & (e > 0) & (f > 0)
        if np.any(mask):
            ax1.loglog(e[mask], e[mask] * f[mask], label=label)
    ax1.set_ylabel(r"$E \times F(E)$")
    ax1.set_title(f"Warm corona slab angle-resolved spectra: {tag}")

    ax2 = fig.add_subplot(2, 1, 2)
    for label, data in angle_data:
        if "poldeg" not in data:
            continue
        e = data["en"]
        p = data["poldeg"]
        mask = np.isfinite(e) & np.isfinite(p) & (e > 0)
        if np.any(mask):
            ax2.semilogx(e[mask], p[mask], label=label)
    ax2.set_xlabel("Energy (keV)")
    ax2.set_ylabel("Polarization degree")

    handles, labels = ax2.get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig(out_dir / f"{tag}_angle_paperstyle.png", dpi=300)
    plt.close()


def make_plots(point_dir: Path, tag: str) -> None:
    """
    读取所有原始 theta 角度 calspec 输出并画图。
    注意：这里保持原功能不变，不读取 mu 分箱结果。
    """
    plot_dir = PLOT_ROOT / tag
    ensure_dir(plot_dir)

    angle_data = []
    for imin in range(0, 90, 3):
        imax = imin + 3
        cdir = point_dir / angle_dir_name(imin, imax)
        if not cdir.exists():
            continue
        if not calspec_outputs_exist(cdir):
            continue

        label = f"{imin:02d}-{imax:02d}"
        data = read_calspec_dir(cdir)
        angle_data.append((label, data))

    if len(angle_data) == 0:
        print("[警告] 没有可用于画图的角度 calspec 结果。")
        return

    plot_angle_flux(angle_data, plot_dir, tag)
    plot_angle_poldeg(angle_data, plot_dir, tag)
    plot_angle_polang(angle_data, plot_dir, tag)
    plot_angle_paperstyle(angle_data, plot_dir, tag)

    print(f"[绘图完成] 输出目录：{plot_dir}")


# ============================================================
# 主流程
# ============================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="自动运行 MONK slab + calspec + angle-resolved polarization plotting")
    parser.add_argument("--te", type=float, required=True, help="electron temperature, keV")
    parser.add_argument("--tau", type=float, required=True, help="optical depth")

    parser.add_argument("--ne", type=int, default=DEFAULT_NE, help="calspec ne, 默认 -10000")
    parser.add_argument("--emin", type=float, default=DEFAULT_EMIN, help="calspec emin, 默认 0.01")
    parser.add_argument("--emax", type=float, default=DEFAULT_EMAX, help="calspec emax, 默认 600.0")

    parser.add_argument("--force-slab", action="store_true", help="强制重跑 slab")
    parser.add_argument("--force-calspec", action="store_true", help="强制重跑总谱 calspec")
    parser.add_argument("--force-angle", action="store_true", help="强制重跑所有 theta 角度 calspec")

    # 新增 mu 功能
    parser.add_argument("--run-mu", action="store_true", help="额外生成 mu=cos(theta) 均匀分箱的 calspec 输出")
    parser.add_argument("--nmu", type=int, default=30, help="mu 分箱个数，默认 30")
    parser.add_argument("--force-mu", action="store_true", help="强制重跑所有 mu 分箱 calspec")

    args = parser.parse_args()

    te = args.te
    tau = args.tau
    ne = args.ne
    emin = args.emin
    emax = args.emax

    if not SLAB_EXE.exists():
        raise FileNotFoundError(f"未找到 slab 程序：{SLAB_EXE}")
    if not CALSPEC_EXE.exists():
        raise FileNotFoundError(f"未找到 calspec 程序：{CALSPEC_EXE}")

    tag = format_tag(te, tau)
    point_dir = BASE_DIR / tag
    slab_dir = point_dir / "slab"

    ensure_dir(point_dir)
    ensure_dir(slab_dir)

    print("=" * 72)
    print("参数点：", tag)
    print("point_dir:", point_dir)
    print("plot_dir :", PLOT_ROOT / tag)
    if args.run_mu:
        print(f"mu_bins  : enabled, nmu={args.nmu}")
    print("=" * 72)

    # 1. 根目录留一份 params.txt 记录
    prepare_params_file(point_dir, te, tau)

    # 2. slab 目录放真正运行使用的 params.txt
    prepare_params_file(slab_dir, te, tau)

    # 3. 跑 slab
    slab_dir = run_slab(point_dir, force=args.force_slab)

    # 4. 复制 q/u -> qweight/uweight
    ensure_qu_as_qweight(slab_dir)

    # 5. 总谱 calspec
    run_calspec_total(
        slab_dir,
        point_dir,
        ne=ne,
        emin=emin,
        emax=emax,
        force=args.force_calspec
    )

    # 6. 原始 theta 角度 calspec（保留原功能）
    for imin in range(0, 90, 3):
        imax = imin + 3
        run_calspec_angle(
            slab_dir,
            point_dir,
            imin=imin,
            imax=imax,
            ne=ne,
            emin=emin,
            emax=emax,
            force=args.force_angle
        )

    # 7. 新增 mu 均匀分箱 calspec（可选）
    if args.run_mu:
        run_all_calspec_mu_bins(
            slab_dir=slab_dir,
            point_dir=point_dir,
            ne=ne,
            emin=emin,
            emax=emax,
            nmu=args.nmu,
            force=args.force_mu,
        )

    # 8. 画图（保持原有 theta 角度绘图逻辑）
    make_plots(point_dir, tag)

    print("=" * 72)
    print("全部完成。")
    print("=" * 72)


if __name__ == "__main__":
    main()