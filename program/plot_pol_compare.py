#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compare_pol_error.py

只针对 MONK slab，做 PA 收敛误差图。

新版逻辑
--------
1. 先生成原始细角度 calspec：
   - theta 模式：00-03, 03-06, ..., 87-90
   - mu 模式：calspec_mu_00 ... calspec_mu_{nmu-1}

2. 后处理时：
   - 按用户要求的大角度区间，把多个原始 calspec 目录的 I/Q/U 直接相加合并
   - 再按 energy edges 对合并后的 I/Q/U 做能量积分
   - 最后从积分后的 Q/U 算 band-integrated PA

3. 对每个 nph，下的 5 次 realization 计算：
   - mean(PA)
   - SE(PA)

4. 对 SE(PA) vs N 拟合：
   SE(PA) = A * N^(-alpha)

输出
----
PLOT_ROOT / te_xxx_tau_xxx /
    pa_convergence.csv
    pa_convergence_single_column.png
    pa_convergence_single_column.pdf
"""

import argparse
import concurrent.futures as cf
import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# 路径设置
# ============================================================
MONK_BIN_DIR = Path("/home/hdw/data/monk/monk_for_rhy/bin")
SLAB_EXE = MONK_BIN_DIR / "slab"
CALSPEC_EXE = MONK_BIN_DIR / "calspec"

PARAM_TEMPLATE = Path("/home/hdw/data/monk/plot/warmcorona/test/polarization/data/10_7/params.txt")

DEFAULT_BASE_ROOT = Path("/home/hdw/data/monk/plot/warmcorona/test/polarization/data/test_slab_err")
DEFAULT_PLOT_ROOT = Path("/home/hdw/data/monk/plot/warmcorona/test/2026.04.02/test_slab_err")

DEFAULT_NE = -10000
DEFAULT_EMIN = 0.01
DEFAULT_EMAX = 600.0


# ============================================================
# 数据结构
# ============================================================
@dataclass
class AngleGroup:
    label: str
    mode: str   # "theta" or "mu"
    v1: float
    v2: float


# ============================================================
# 基础工具
# ============================================================
def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def run_bash(command: str, cwd: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(
        ["bash", "-lc", command],
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True
    )
    if check and result.returncode != 0:
        print("\n[ERROR] shell command failed")
        print(command)
        print("stdout:\n", result.stdout)
        print("stderr:\n", result.stderr)
        raise RuntimeError("shell command failed")
    return result


def save_log(path: Path, title: str, result: subprocess.CompletedProcess) -> None:
    text = []
    text.append(f"COMMAND:\n{title}\n")
    text.append("STDOUT:\n")
    text.append(result.stdout)
    text.append("\nSTDERR:\n")
    text.append(result.stderr)
    path.write_text("\n".join(text), encoding="utf-8")


def format_tag(te: float, tau: float) -> str:
    return f"te_{te:.3f}_tau_{tau:.3f}"


def photon_label(nph: int) -> str:
    """
    100000 -> 10_5
    """
    if nph <= 0:
        raise ValueError("nph must be positive")
    exp = int(round(np.log10(float(nph))))
    if not np.isclose(nph, 10 ** exp):
        return str(nph)
    return f"10_{exp}"


def load_dat_file(file_path: Path) -> np.ndarray:
    raw = file_path.read_bytes()
    if len(raw) % 8 != 0:
        raise ValueError(f"{file_path} 不是 8 字节倍数，无法按 little-endian double 解析")
    return np.frombuffer(raw, dtype="<f8").copy()


def geometric_centers_from_edges(edges: np.ndarray) -> np.ndarray:
    return np.sqrt(edges[:-1] * edges[1:])


# ============================================================
# params.txt 处理
# ============================================================
def prepare_params_file(
    out_path: Path,
    te: float,
    tau: float,
    nph: int,
    template_path: Path,
    photon_keys: List[str]
) -> None:
    if not template_path.exists():
        raise FileNotFoundError(f"未找到 params 模板: {template_path}")

    text = template_path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    te_done = False
    tau_done = False
    photon_done = False
    new_lines = []

    for line in lines:
        stripped = line.strip()

        if re.match(r"^te\b", stripped):
            new_lines.append(f"te = {te:.3f}")
            te_done = True
            continue

        if re.match(r"^tau\b", stripped):
            new_lines.append(f"tau = {tau:.3f}")
            tau_done = True
            continue

        matched_photon = False
        for key in photon_keys:
            if re.match(rf"^{re.escape(key)}\b", stripped):
                new_lines.append(f"{key} = {nph}")
                photon_done = True
                matched_photon = True
                break
        if matched_photon:
            continue

        new_lines.append(line)

    if not te_done:
        new_lines.append(f"te = {te:.3f}")
    if not tau_done:
        new_lines.append(f"tau = {tau:.3f}")
    if not photon_done:
        new_lines.append(f"{photon_keys[0]} = {nph}")

    out_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")


# ============================================================
# slab / calspec 检查
# ============================================================
def slab_outputs_exist(slab_dir: Path) -> bool:
    return (slab_dir / "en0.dat").exists() and (slab_dir / "weight.dat").exists()


def calspec_outputs_exist(cdir: Path) -> bool:
    return (cdir / "en.dat").exists() and (cdir / "de.dat").exists() and (cdir / "flux.dat").exists()


def ensure_qweight_uweight(slab_dir: Path) -> None:
    qarr = slab_dir / "qarr.dat"
    uarr = slab_dir / "uarr.dat"
    qweight = slab_dir / "qweight.dat"
    uweight = slab_dir / "uweight.dat"

    if qarr.exists() and not qweight.exists():
        shutil.copy2(qarr, qweight)
    if uarr.exists() and not uweight.exists():
        shutil.copy2(uarr, uweight)


# ============================================================
# 角度 group 解析（用户输入的大区间）
# ============================================================
def parse_theta_groups(spec: str) -> List[AngleGroup]:
    groups = []
    items = [x.strip() for x in spec.split(",") if x.strip()]
    for item in items:
        a, b = item.split("-")
        amin = float(a)
        amax = float(b)
        label = f"theta_{int(round(amin)):02d}_{int(round(amax)):02d}"
        groups.append(AngleGroup(label=label, mode="theta", v1=amin, v2=amax))
    return groups


def parse_mu_groups(spec: str) -> List[AngleGroup]:
    groups = []
    items = [x.strip() for x in spec.split(",") if x.strip()]
    for item in items:
        a, b = item.split("-")
        mu1 = float(a)
        mu2 = float(b)
        label = f"mu_{mu1:.3f}_{mu2:.3f}".replace(".", "p")
        groups.append(AngleGroup(label=label, mode="mu", v1=mu1, v2=mu2))
    return groups


# ============================================================
# 原始 native bins
# ============================================================
def native_theta_bins() -> List[Tuple[int, int]]:
    out = []
    for imin in range(0, 90, 3):
        imax = imin + 3
        out.append((imin, imax))
    return out


def native_theta_dir_name(imin: int, imax: int) -> str:
    return f"calspec_{imin:02d}-{imax:02d}"


def native_mu_dir_name(idx: int) -> str:
    return f"calspec_mu_{idx:02d}"


def theta_group_to_native_dirs(run_dir: Path, theta_min: float, theta_max: float) -> List[Path]:
    """
    把大角度区间映射到原始 3° bins，并返回对应目录列表。
    例如 0-12 -> 00-03,03-06,06-09,09-12
    """
    dirs = []
    for imin, imax in native_theta_bins():
        if imin >= theta_min and imax <= theta_max:
            cdir = run_dir / native_theta_dir_name(imin, imax)
            if cdir.exists() and calspec_outputs_exist(cdir):
                dirs.append(cdir)
    return dirs


def mu_group_to_native_dirs(run_dir: Path, mu_min: float, mu_max: float, nmu: int) -> List[Path]:
    """
    把大 mu 区间映射到原始 mu bins。
    原始 mu bin 定义为 [j/nmu, (j+1)/nmu]
    这里取完全落在 [mu_min, mu_max] 内的 bins。
    """
    dirs = []
    edges = np.linspace(0.0, 1.0, nmu + 1)
    for j in range(nmu):
        left = float(edges[j])
        right = float(edges[j + 1])
        if left >= mu_min and right <= mu_max:
            cdir = run_dir / native_mu_dir_name(j)
            if cdir.exists() and calspec_outputs_exist(cdir):
                dirs.append(cdir)
    return dirs


# ============================================================
# 能量区间
# ============================================================
def build_energy_edges(
    e_min: float,
    e_max: float,
    nbins: int,
    explicit_edges: Optional[str] = None
) -> np.ndarray:
    if explicit_edges:
        arr = np.array([float(x.strip()) for x in explicit_edges.split(",") if x.strip()], dtype=float)
        if len(arr) < 2:
            raise ValueError("显式 energy edges 至少要有两个值")
        if not np.all(arr[1:] > arr[:-1]):
            raise ValueError("energy edges 必须严格递增")
        return arr

    if e_min <= 0 or e_max <= e_min or nbins <= 0:
        raise ValueError("非法能量分箱参数")
    return np.logspace(np.log10(e_min), np.log10(e_max), nbins + 1)


# ============================================================
# 偏振角 / 统计
# ============================================================
def normalize_pa_deg(pa_deg: float) -> float:
    return ((pa_deg + 90.0) % 180.0) - 90.0


def pa_from_stokes(Q: float, U: float) -> float:
    pa = 0.5 * math.degrees(math.atan2(U, Q))
    return normalize_pa_deg(pa)


def circular_mean_and_se_pa(pa_list_deg: List[float]) -> Tuple[float, float]:
    vals = np.array([x for x in pa_list_deg if np.isfinite(x)], dtype=float)
    n = len(vals)
    if n == 0:
        return np.nan, np.nan
    if n == 1:
        return normalize_pa_deg(float(vals[0])), np.nan

    phi = np.deg2rad(2.0 * vals)
    c = np.mean(np.cos(phi))
    s = np.mean(np.sin(phi))
    mean_phi = math.atan2(s, c)
    mean_pa = normalize_pa_deg(0.5 * math.degrees(mean_phi))

    mean_phi_arr = np.deg2rad(2.0 * np.full(n, mean_pa))
    dphi = phi - mean_phi_arr
    dphi = (dphi + np.pi) % (2.0 * np.pi) - np.pi
    dpa_deg = np.rad2deg(dphi) / 2.0

    se = np.sqrt(np.sum(dpa_deg ** 2) / (n * (n - 1)))
    return mean_pa, float(se)


def fit_power_law_alpha(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.sum(mask) < 2:
        return np.nan, np.nan

    lx = np.log10(x[mask])
    ly = np.log10(y[mask])

    coeff = np.polyfit(lx, ly, 1)
    slope, intercept = coeff[0], coeff[1]
    alpha = -slope
    A = 10.0 ** intercept
    return float(alpha), float(A)


# ============================================================
# Stokes 积分 / 合并
# ============================================================
def integrate_stokes_band(
    E: np.ndarray,
    dE: np.ndarray,
    I: np.ndarray,
    Q: np.ndarray,
    U: np.ndarray,
    emin: float,
    emax: float
) -> Tuple[float, float, float, float]:
    """
    返回 I_int, Q_int, U_int, PA_int
    """
    mask = (
        np.isfinite(E) & np.isfinite(dE) &
        np.isfinite(I) & np.isfinite(Q) & np.isfinite(U) &
        (dE > 0.0) &
        (E >= emin) & (E <= emax)
    )

    if not np.any(mask):
        return np.nan, np.nan, np.nan, np.nan

    I_int = np.sum(I[mask] * dE[mask])
    Q_int = np.sum(Q[mask] * dE[mask])
    U_int = np.sum(U[mask] * dE[mask])

    if not np.isfinite(I_int) or I_int <= 0:
        return I_int, Q_int, U_int, np.nan

    pa = pa_from_stokes(Q_int, U_int)
    return float(I_int), float(Q_int), float(U_int), float(pa)


def integrate_stokes_over_edges(
    E: np.ndarray,
    dE: np.ndarray,
    I: np.ndarray,
    Q: np.ndarray,
    U: np.ndarray,
    edges: np.ndarray
) -> List[float]:
    pa_list = []
    for i in range(len(edges) - 1):
        left = float(edges[i])
        right = float(edges[i + 1])
        _, _, _, pa = integrate_stokes_band(E, dE, I, Q, U, left, right)
        pa_list.append(pa)
    return pa_list


def merge_multiple_calspec_dirs(cdirs: List[Path]) -> Optional[Dict[str, np.ndarray]]:
    """
    像 plot_pol_compare.py 那样：
    直接把多个 calspec 目录在同一能量网格上的 I/Q/U 相加
    """
    valid_dirs = []
    for c in cdirs:
        qf = c / "qflux.dat"
        uf = c / "uflux.dat"
        if c.exists() and calspec_outputs_exist(c) and qf.exists() and uf.exists():
            valid_dirs.append(c)

    if len(valid_dirs) == 0:
        return None

    E0 = load_dat_file(valid_dirs[0] / "en.dat")
    dE0 = load_dat_file(valid_dirs[0] / "de.dat")
    I_sum = np.zeros_like(E0)
    Q_sum = np.zeros_like(E0)
    U_sum = np.zeros_like(E0)

    for cdir in valid_dirs:
        E = load_dat_file(cdir / "en.dat")
        dE = load_dat_file(cdir / "de.dat")
        I = load_dat_file(cdir / "flux.dat")
        Q = load_dat_file(cdir / "qflux.dat")
        U = load_dat_file(cdir / "uflux.dat")

        if len(E) != len(E0) or not np.allclose(E, E0, rtol=1e-10, atol=1e-12):
            raise RuntimeError(f"能量网格不一致，无法合并：{cdir}")
        if len(dE) != len(dE0) or not np.allclose(dE, dE0, rtol=1e-10, atol=1e-12):
            raise RuntimeError(f"dE 网格不一致，无法合并：{cdir}")

        I_sum += I
        Q_sum += Q
        U_sum += U

    return {
        "E": E0.copy(),
        "dE": dE0.copy(),
        "I": I_sum,
        "Q": Q_sum,
        "U": U_sum,
    }


# ============================================================
# 运行一个 realization：只生成原始 native bins
# ============================================================
def run_one_realization(
    te: float,
    tau: float,
    nph: int,
    run_idx: int,
    point_dir: Path,
    angle_mode: str,
    nmu: int,
    ne: int,
    emin: float,
    emax: float,
    template_path: Path,
    force_slab: bool,
    force_calspec: bool,
    photon_keys: List[str]
) -> Dict[str, str]:
    run_dir = point_dir / f"run_{run_idx:02d}"
    slab_dir = run_dir / "slab"
    ensure_dir(slab_dir)

    prepare_params_file(
        out_path=slab_dir / "params.txt",
        te=te,
        tau=tau,
        nph=nph,
        template_path=template_path,
        photon_keys=photon_keys
    )

    if (not slab_outputs_exist(slab_dir)) or force_slab:
        if not SLAB_EXE.exists():
            raise FileNotFoundError(f"未找到 slab 程序: {SLAB_EXE}")
        cmd = f"cd '{slab_dir}' && '{SLAB_EXE}'"
        result = run_bash(cmd, check=True)
        save_log(slab_dir / "slab.log", cmd, result)

    ensure_qweight_uweight(slab_dir)

    if not CALSPEC_EXE.exists():
        raise FileNotFoundError(f"未找到 calspec 程序: {CALSPEC_EXE}")

    if angle_mode == "theta":
        for imin, imax in native_theta_bins():
            cdir = run_dir / native_theta_dir_name(imin, imax)
            ensure_dir(cdir)

            if calspec_outputs_exist(cdir) and (not force_calspec):
                continue

            cmd = (
                f"cd '{cdir}' && "
                f"'{CALSPEC_EXE}' '{slab_dir}' {ne} {emin} {emax} {imin} {imax}"
            )
            result = run_bash(cmd, check=True)
            save_log(cdir / "calspec.log", cmd, result)

    elif angle_mode == "mu":
        mu_edges = np.linspace(0.0, 1.0, nmu + 1)
        for j in range(nmu):
            mu_min = float(mu_edges[j])
            mu_max = float(mu_edges[j + 1])
            theta_min = math.degrees(math.acos(mu_max))
            theta_max = math.degrees(math.acos(mu_min))

            cdir = run_dir / native_mu_dir_name(j)
            ensure_dir(cdir)

            if calspec_outputs_exist(cdir) and (not force_calspec):
                continue

            cmd = (
                f"cd '{cdir}' && "
                f"'{CALSPEC_EXE}' '{slab_dir}' {ne} {emin} {emax} {theta_min:.8f} {theta_max:.8f}"
            )
            result = run_bash(cmd, check=True)
            save_log(cdir / "calspec.log", cmd, result)

    else:
        raise ValueError(f"未知 angle_mode: {angle_mode}")

    return {
        "run_dir": str(run_dir),
        "slab_dir": str(slab_dir),
        "status": "ok"
    }


def worker_run_one_realization(payload):
    return run_one_realization(
        te=payload["te"],
        tau=payload["tau"],
        nph=payload["nph"],
        run_idx=payload["run_idx"],
        point_dir=Path(payload["point_dir"]),
        angle_mode=payload["angle_mode"],
        nmu=payload["nmu"],
        ne=payload["ne"],
        emin=payload["emin"],
        emax=payload["emax"],
        template_path=Path(payload["template_path"]),
        force_slab=payload["force_slab"],
        force_calspec=payload["force_calspec"],
        photon_keys=payload["photon_keys"],
    )


# ============================================================
# 从一个 realization 中读取“大角区间 × 能量区间”的 PA
# ============================================================
def read_pa_for_group_from_run(
    run_dir: Path,
    group: AngleGroup,
    nmu: int,
    edges: np.ndarray
) -> List[float]:
    if group.mode == "theta":
        cdirs = theta_group_to_native_dirs(run_dir, group.v1, group.v2)
    elif group.mode == "mu":
        cdirs = mu_group_to_native_dirs(run_dir, group.v1, group.v2, nmu)
    else:
        raise ValueError(f"未知 group mode: {group.mode}")

    merged = merge_multiple_calspec_dirs(cdirs)
    if merged is None:
        return [np.nan] * (len(edges) - 1)

    return integrate_stokes_over_edges(
        merged["E"], merged["dE"], merged["I"], merged["Q"], merged["U"], edges
    )


# ============================================================
# 绘图
# ============================================================
def plot_single_column(
    df: pd.DataFrame,
    groups: List[AngleGroup],
    edges: np.ndarray,
    out_png: Path,
    out_pdf: Path,
    alpha_threshold: float = 0.35,
    pa_err_ref: Optional[float] = 4.0
) -> None:
    nrow = len(groups)
    fig, axes = plt.subplots(nrow, 1, figsize=(8.2, 3.4 * nrow), sharex=True)
    if nrow == 1:
        axes = [axes]

    centers = geometric_centers_from_edges(edges)
    color_cycle = ["green", "black", "red", "blue", "purple", "orange", "brown", "deepskyblue"]

    for ia, group in enumerate(groups):
        ax = axes[ia]
        sub = df[df["angle_label"] == group.label].copy()

        has_positive_y = False
        has_any_curve = False

        for ibin in range(len(edges) - 1):
            s2 = sub[sub["ebin"] == ibin].sort_values("nph")
            if len(s2) == 0:
                continue

            x = s2["nph"].values.astype(float)
            y = s2["pa_se_deg"].values.astype(float)

            mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            if np.sum(mask) == 0:
                continue

            x_plot = x[mask]
            y_plot = y[mask]

            has_any_curve = True
            has_positive_y = True

            alpha = s2["alpha_fit"].iloc[0] if len(s2) > 0 else np.nan
            afit = s2["A_fit"].iloc[0] if len(s2) > 0 else np.nan

            color = color_cycle[ibin % len(color_cycle)]
            label = f"{centers[ibin]:.2f} keV  {alpha:.2f}"

            ax.plot(x_plot, y_plot, marker="o", color=color, lw=1.6, label=label)

            if np.isfinite(alpha) and np.isfinite(afit) and len(x_plot) >= 2:
                xx = np.logspace(np.log10(np.min(x_plot)), np.log10(np.max(x_plot)), 100)
                yy = afit * xx ** (-alpha)
                good = np.isfinite(yy) & (yy > 0)
                if np.any(good):
                    ax.plot(xx[good], yy[good], color=color, alpha=0.45, lw=2.5)

        if pa_err_ref is not None and pa_err_ref > 0:
            ax.axhline(pa_err_ref, color="gray", lw=1.0, ls="--", alpha=0.9)

        if group.mode == "theta":
            ax.set_title(f"{group.v1:.0f}–{group.v2:.0f} deg", fontsize=13)
        else:
            ax.set_title(f"mu={group.v1:.2f}–{group.v2:.2f}", fontsize=13)

        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(matplotlib.ticker.LogFormatterMathtext())

        if has_positive_y:
            ax.set_yscale("log")
        else:
            ax.set_yscale("linear")
            ax.set_ylim(0, 1)
            ax.text(
                0.5, 0.5, "No positive SE(PA) values",
                transform=ax.transAxes,
                ha="center", va="center",
                fontsize=11,
                color="gray"
            )

        ax.set_ylabel("SE of PA (deg)", fontsize=12)
        ax.grid(alpha=0.25, which="both")

        if has_any_curve:
            ax.legend(fontsize=10, frameon=True, loc="best")

        bad = sub[["ebin", "alpha_fit"]].drop_duplicates()
        bad_list = []
        for _, row in bad.iterrows():
            if np.isfinite(row["alpha_fit"]) and row["alpha_fit"] < alpha_threshold:
                ebin = int(row["ebin"])
                bad_list.append(f"{centers[ebin]:.2f}keV")

        if bad_list:
            txt = r"$\alpha < $" + f"{alpha_threshold:.2f}: " + ", ".join(bad_list)
            ax.text(
                0.98, 0.05, txt,
                ha="right", va="bottom",
                transform=ax.transAxes,
                fontsize=10,
                bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.75)
            )

    axes[-1].set_xlabel("Super-photon number", fontsize=13)
    fig.subplots_adjust(left=0.12, right=0.96, top=0.96, bottom=0.08, hspace=0.35)

    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


# ============================================================
# 主流程
# ============================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="slab PA convergence test with native-angle merging")

    parser.add_argument("--te", type=float, required=True, help="electron temperature")
    parser.add_argument("--tau", type=float, required=True, help="optical depth")

    parser.add_argument("--base-root", type=str, default=str(DEFAULT_BASE_ROOT), help="run output root")
    parser.add_argument("--plot-root", type=str, default=str(DEFAULT_PLOT_ROOT), help="plot/csv output root")
    parser.add_argument("--template", type=str, default=str(PARAM_TEMPLATE), help="params.txt template")

    parser.add_argument("--nphotons", type=str, default="1e5,1e6,1e7",
                        help="comma-separated photon numbers, e.g. 1e5,1e6,1e7")
    parser.add_argument("--nruns", type=int, default=5, help="independent runs per photon number")
    parser.add_argument("--jobs", type=int, default=3, help="parallel workers")

    parser.add_argument("--ne", type=int, default=DEFAULT_NE)
    parser.add_argument("--emin", type=float, default=DEFAULT_EMIN)
    parser.add_argument("--emax", type=float, default=DEFAULT_EMAX)

    parser.add_argument("--angle-mode", choices=["theta", "mu"], default="theta")
    parser.add_argument("--theta-bins", type=str, default="0-12,30-40,60-70,80-90",
                        help="when angle-mode=theta, e.g. 0-12,30-40,60-70,80-90")
    parser.add_argument("--mu-bins", type=str, default="0.0-0.2,0.2-0.5,0.5-0.8,0.8-1.0",
                        help="when angle-mode=mu, e.g. 0.0-0.2,0.2-0.5,0.5-0.8,0.8-1.0")
    parser.add_argument("--nmu", type=int, default=30, help="number of native mu bins")

    parser.add_argument("--band-emin", type=float, default=1.5, help="energy-integration range min")
    parser.add_argument("--band-emax", type=float, default=10.0, help="energy-integration range max")
    parser.add_argument("--band-nbins", type=int, default=3, help="number of integrated energy bins")
    parser.add_argument("--energy-edges", type=str, default="",
                        help='explicit energy edges, e.g. "1.5,2.5,5.0,10.0"')

    parser.add_argument("--force-slab", action="store_true", help="force rerun slab")
    parser.add_argument("--force-calspec", action="store_true", help="force rerun calspec")
    parser.add_argument("--alpha-threshold", type=float, default=0.35,
                        help="alpha threshold for PA reliability")
    parser.add_argument("--pa-err-ref", type=float, default=4.0,
                        help="reference horizontal line in SE(PA) plot")

    parser.add_argument("--photon-keys", type=str, default="nph,nphoton,nsuperphoton",
                        help="comma-separated possible parameter names in params.txt")

    args = parser.parse_args()

    te = args.te
    tau = args.tau
    base_root = Path(args.base_root)
    plot_root = Path(args.plot_root)
    template_path = Path(args.template)
    photon_keys = [x.strip() for x in args.photon_keys.split(",") if x.strip()]

    if args.angle_mode == "theta":
        groups = parse_theta_groups(args.theta_bins)
    else:
        groups = parse_mu_groups(args.mu_bins)

    edges = build_energy_edges(
        e_min=args.band_emin,
        e_max=args.band_emax,
        nbins=args.band_nbins,
        explicit_edges=args.energy_edges if args.energy_edges.strip() else None
    )

    nph_list = []
    for s in [x.strip() for x in args.nphotons.split(",") if x.strip()]:
        nph_list.append(int(float(s)))
    nph_list = sorted(set(nph_list))

    tag = format_tag(te, tau)
    plot_dir = plot_root / tag
    ensure_dir(plot_dir)

    tasks = []
    for nph in nph_list:
        point_dir = base_root / photon_label(nph) / tag
        ensure_dir(point_dir)
        for irun in range(1, args.nruns + 1):
            tasks.append({
                "te": te,
                "tau": tau,
                "nph": nph,
                "run_idx": irun,
                "point_dir": str(point_dir),
                "angle_mode": args.angle_mode,
                "nmu": args.nmu,
                "ne": args.ne,
                "emin": args.emin,
                "emax": args.emax,
                "template_path": str(template_path),
                "force_slab": args.force_slab,
                "force_calspec": args.force_calspec,
                "photon_keys": photon_keys,
            })

    print("=" * 72)
    print("[INFO] tag       :", tag)
    print("[INFO] base_root :", base_root)
    print("[INFO] plot_root :", plot_root)
    print("[INFO] nph list   :", nph_list)
    print("[INFO] nruns      :", args.nruns)
    print("[INFO] jobs       :", args.jobs)
    print("[INFO] angle mode :", args.angle_mode)
    print("[INFO] groups     :", [g.label for g in groups])
    print("[INFO] E edges    :", edges.tolist())
    print("=" * 72)

    if args.jobs <= 1:
        for t in tasks:
            worker_run_one_realization(t)
    else:
        with cf.ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futures = [ex.submit(worker_run_one_realization, t) for t in tasks]
            done = 0
            total = len(futures)
            for fu in cf.as_completed(futures):
                fu.result()
                done += 1
                print(f"[DONE {done}/{total}]")

    records = []
    for nph in nph_list:
        point_dir = base_root / photon_label(nph) / tag

        for group in groups:
            pa_runs: Dict[int, List[float]] = {i: [] for i in range(len(edges) - 1)}

            for irun in range(1, args.nruns + 1):
                run_dir = point_dir / f"run_{irun:02d}"

                try:
                    pa_list = read_pa_for_group_from_run(run_dir, group, args.nmu, edges)
                    for i, pa in enumerate(pa_list):
                        pa_runs[i].append(pa)
                except Exception as e:
                    print(f"[WARN] failed to read {run_dir} / {group.label}: {e}")

            for ebin in range(len(edges) - 1):
                mean_pa, se_pa = circular_mean_and_se_pa(pa_runs[ebin])
                records.append({
                    "tag": tag,
                    "te": te,
                    "tau": tau,
                    "angle_mode": group.mode,
                    "angle_label": group.label,
                    "angle_v1": group.v1,
                    "angle_v2": group.v2,
                    "ebin": ebin,
                    "e1": float(edges[ebin]),
                    "e2": float(edges[ebin + 1]),
                    "ecenter": float(math.sqrt(edges[ebin] * edges[ebin + 1])),
                    "nph": nph,
                    "nruns_used": len([x for x in pa_runs[ebin] if np.isfinite(x)]),
                    "pa_mean_deg": mean_pa,
                    "pa_se_deg": se_pa,
                })

    df = pd.DataFrame.from_records(records)
    if len(df) == 0:
        raise RuntimeError("没有任何统计结果，检查运行目录是否正确。")

    alpha_records = []
    group_cols = ["angle_label", "ebin"]
    for keys, sub in df.groupby(group_cols):
        x = sub["nph"].values.astype(float)
        y = sub["pa_se_deg"].values.astype(float)
        alpha, A = fit_power_law_alpha(x, y)
        alpha_records.append({
            "angle_label": keys[0],
            "ebin": keys[1],
            "alpha_fit": alpha,
            "A_fit": A
        })

    df_alpha = pd.DataFrame(alpha_records)
    df = df.merge(df_alpha, on=["angle_label", "ebin"], how="left")

    csv_path = plot_dir / "pa_convergence.csv"
    df.sort_values(["angle_label", "ebin", "nph"]).to_csv(csv_path, index=False)

    out_png = plot_dir / "pa_convergence_single_column.png"
    out_pdf = plot_dir / "pa_convergence_single_column.pdf"

    plot_single_column(
        df=df,
        groups=groups,
        edges=edges,
        out_png=out_png,
        out_pdf=out_pdf,
        alpha_threshold=args.alpha_threshold,
        pa_err_ref=args.pa_err_ref
    )

    print("=" * 72)
    print("[OK] CSV :", csv_path)
    print("[OK] PNG :", out_png)
    print("[OK] PDF :", out_pdf)
    print("=" * 72)


if __name__ == "__main__":
    main()