#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import math
import random
import subprocess
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =========================================================
# 1. 路径
# =========================================================

BASE_EXPORT_DIR = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"

OUT_ROOT = os.path.expanduser(
    "~/data/monk/plot/warmcorona/test/2026.03.25/"
)
os.makedirs(OUT_ROOT, exist_ok=True)

# -------- local model ----------
LOCAL_MODEL_DIR  = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_slab_log/warmcom_slab"
LOCAL_PKG_NAME   = "warmcom_slab"
LOCAL_MODEL_NAME = "warmcomslab"

# -------- mod model ----------
MOD_MODEL_PATH = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_pure_Tom_slab_107.mod"


# =========================================================
# 2. 参数
# =========================================================

# ---- 节点对比：这些点必须在 _exports 目录真实存在 ----
NODE_PARAM_SETS = [
    (0.300, 9.500),
    (0.300, 15.000),
    (0.500, 19.500),
    (0.600, 9.800),
]

# ---- 随机插值点对比：不要求在 _exports 中存在 ----
RANDOM_COMPARE_N = 8
RANDOM_SEED = 12345

# 参数空间范围
TE_MIN = 0.1
TE_MAX = 0.6
TAU_MIN = 5.0
TAU_MAX = 25.0

# 是否避开整数/明显节点，专门测“插值点”
AVOID_GRIDLIKE_POINTS = True

# 模型公共参数
Z = 0.0
NORM = 1.0

# 归一化参考能量
E_REF = 0.15

# XSPEC plot energy grid
E_MIN = 0.001
E_MAX = 500.0
NBINS = 1000
GRID_MODE = "log"

SAVE_XCM = True
PRINT_DEBUG = True


# =========================================================
# 3. 工具函数
# =========================================================

def read_binary_double(path: str) -> np.ndarray:
    arr = np.fromfile(path, dtype="<d")
    if arr.size == 0:
        raise RuntimeError(f"读取失败或空文件: {path}")
    return arr


def read_qdp_model(path: str):
    xs, ys = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if (not s) or s[0] in ("!", "@"):
                continue

            parts = s.split()
            try:
                if len(parts) >= 3:
                    x = float(parts[0])
                    y = float(parts[2])
                elif len(parts) >= 2:
                    x = float(parts[0])
                    y = float(parts[1])
                else:
                    continue
            except Exception:
                continue

            xs.append(x)
            ys.append(y)

    if not xs:
        raise RuntimeError(f"QDP 为空: {path}")

    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[ok]
    y = y[ok]

    if x.size < 2:
        raise RuntimeError(f"QDP 有效点不足: {path}")

    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    return x, y


def interp_loglog(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.count_nonzero(ok) < 2:
        raise RuntimeError("log-log 插值失败：有效点不足")

    x_ok = x[ok]
    y_ok = y[ok]

    idx = np.argsort(x_ok)
    x_ok = x_ok[idx]
    y_ok = y_ok[idx]

    lx = np.log10(x_ok)
    ly = np.log10(y_ok)
    lx_new = np.log10(x_new)

    ly_new = np.interp(lx_new, lx, ly)
    return 10.0 ** ly_new


def normalize_at(E: np.ndarray, F: np.ndarray, e_ref: float) -> np.ndarray:
    ok = np.isfinite(E) & np.isfinite(F) & (E > 0) & (F > 0)
    if np.count_nonzero(ok) < 2:
        raise RuntimeError("归一化失败：有效点不足")

    f_ref = interp_loglog(np.array([e_ref]), E[ok], F[ok])[0]
    if (not np.isfinite(f_ref)) or (f_ref <= 0):
        raise RuntimeError(f"在 {e_ref} keV 处归一化失败")

    return F / f_ref


def run_xspec_script(script_text: str, workdir: str, extra_env=None):
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env.setdefault("PGPLOT_DEV", "/null")

    if extra_env is not None:
        env.update(extra_env)

    proc = subprocess.run(
        ["xspec"],
        input=script_text,
        text=True,
        cwd=workdir,
        env=env,
        capture_output=True
    )

    if proc.returncode != 0:
        raise RuntimeError(
            "XSPEC 运行失败\n"
            "---------------- STDOUT tail ----------------\n"
            f"{proc.stdout[-5000:]}\n"
            "---------------- STDERR tail ----------------\n"
            f"{proc.stderr[-5000:]}\n"
        )
    return proc


def export_local_model_qdp(te: float, tau: float, out_qdp: str, out_xcm: str | None = None):
    script = f"""
query yes
cpd /null
lmod {LOCAL_PKG_NAME} {LOCAL_MODEL_DIR}
model {LOCAL_MODEL_NAME}
{te}
{tau}
{Z}
{NORM}
setplot rebin 1 1
setplot energy
setplot area off
energies {E_MIN} {E_MAX} {NBINS} {GRID_MODE}
plot model
iplot
wdata {os.path.basename(out_qdp)}
quit
exit
""".strip() + "\n"

    if out_xcm is not None:
        with open(out_xcm, "w", encoding="utf-8") as f:
            f.write(script)

    run_xspec_script(
        script,
        os.path.dirname(out_qdp),
        extra_env={"WARMCOM_BASEDIR": BASE_EXPORT_DIR}
    )

    if (not os.path.exists(out_qdp)) or os.path.getsize(out_qdp) == 0:
        raise RuntimeError(f"local QDP 未生成: {out_qdp}")


def export_mod_model_qdp(te: float, tau: float, out_qdp: str, out_xcm: str | None = None):
    script = f"""
query yes
cpd /null
model atable{{{MOD_MODEL_PATH}}}
{te}
{tau}
{Z}
{NORM}
setplot rebin 1 1
setplot energy
setplot area off
energies {E_MIN} {E_MAX} {NBINS} {GRID_MODE}
plot model
iplot
wdata {os.path.basename(out_qdp)}
quit
exit
""".strip() + "\n"

    if out_xcm is not None:
        with open(out_xcm, "w", encoding="utf-8") as f:
            f.write(script)

    run_xspec_script(script, os.path.dirname(out_qdp), extra_env=None)

    if (not os.path.exists(out_qdp)) or os.path.getsize(out_qdp) == 0:
        raise RuntimeError(f"mod QDP 未生成: {out_qdp}")


def get_calspec_dir(te: float, tau: float) -> str:
    dirname = f"te_{te:.3f}_tau_{tau:.3f}"
    return os.path.join(BASE_EXPORT_DIR, dirname, "calspec")


def print_debug_block(name: str, E: np.ndarray, F: np.ndarray):
    print(f"  [{name}]")
    print(f"    N = {E.size}")
    print(f"    E min/max = {np.nanmin(E):.6g} / {np.nanmax(E):.6g}")
    print(f"    F min/max = {np.nanmin(F):.6g} / {np.nanmax(F):.6g}")


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def make_node_plot(
    E_raw, F_raw_n, F_local_n, F_mod_n,
    ratio_local, ratio_mod,
    te, tau, out_png, out_pdf
):
    fig = plt.figure(figsize=(8, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.2], hspace=0.05)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

    ax1.loglog(
        E_raw, F_raw_n,
        lw=0,
        marker="o",
        ms=3.0,
        color="black",
        alpha=0.9,
        label="flux_smoothed_tv"
    )
    ax1.loglog(
        E_raw, F_local_n,
        lw=2.2,
        color="red",
        linestyle="-",
        label="Local model"
    )
    ax1.loglog(
        E_raw, F_mod_n,
        lw=2.2,
        color="blue",
        linestyle="--",
        label="Mod model"
    )

    ax1.set_ylabel("Normalized flux")
    ax1.set_title(f"Node compare: te={te:.3f}, tau={tau:.3f}, normalized at {E_REF} keV")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(frameon=False)

    ax2.semilogx(E_raw, ratio_local, lw=1.8, color="red", label="Local / TV")
    ax2.semilogx(E_raw, ratio_mod,   lw=1.8, color="blue", linestyle="--", label="Mod / TV")
    ax2.axhline(1.0, color="black", ls="--", lw=1)

    ax2.set_xlabel("Energy [keV]")
    ax2.set_ylabel("Ratio")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(frameon=False)

    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def make_random_plot(
    E, F_local_n, F_mod_n, ratio_mod_over_local,
    te, tau, out_png, out_pdf
):
    fig = plt.figure(figsize=(8, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.2], hspace=0.05)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

    ax1.loglog(
        E, F_local_n,
        lw=2.2,
        color="red",
        linestyle="-",
        label="Local model"
    )
    ax1.loglog(
        E, F_mod_n,
        lw=2.2,
        color="blue",
        linestyle="--",
        label="Mod model"
    )

    ax1.set_ylabel("Normalized flux")
    ax1.set_title(f"Random interpolation compare: te={te:.4f}, tau={tau:.4f}, normalized at {E_REF} keV")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(frameon=False)

    ax2.semilogx(E, ratio_mod_over_local, lw=1.8, color="purple", label="Mod / Local")
    ax2.axhline(1.0, color="black", ls="--", lw=1)

    ax2.set_xlabel("Energy [keV]")
    ax2.set_ylabel("Ratio")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(frameon=False)

    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def discover_existing_nodes(base_dir: str):
    nodes = []
    pat = re.compile(r"^te_([0-9.]+)_tau_([0-9.]+)$")
    for name in os.listdir(base_dir):
        m = pat.match(name)
        if not m:
            continue
        te = float(m.group(1))
        tau = float(m.group(2))
        calspec_dir = os.path.join(base_dir, name, "calspec")
        if os.path.isdir(calspec_dir):
            nodes.append((te, tau))
    nodes = sorted(set(nodes))
    return nodes


def is_close_to_any_grid_value(x: float, grid_values, tol=1e-8):
    for g in grid_values:
        if abs(x - g) < tol:
            return True
    return False


def choose_random_param_sets(existing_nodes, n, seed=12345, avoid_gridlike=True):
    rng = random.Random(seed)

    te_grid = sorted(set(te for te, _ in existing_nodes))
    tau_grid = sorted(set(tau for _, tau in existing_nodes))

    picked = []
    used = set()

    attempts = 0
    max_attempts = 100000

    while len(picked) < n and attempts < max_attempts:
        attempts += 1

        te = rng.uniform(TE_MIN, TE_MAX)
        tau = rng.uniform(TAU_MIN, TAU_MAX)

        # 保留三位小数，和你的目录/参数风格一致
        te = round(te, 3)
        tau = round(tau, 3)

        if avoid_gridlike:
            if is_close_to_any_grid_value(te, te_grid, tol=1e-12):
                continue
            if is_close_to_any_grid_value(tau, tau_grid, tol=1e-12):
                continue

        key = (te, tau)
        if key in used:
            continue

        used.add(key)
        picked.append(key)

    if len(picked) < n:
        raise RuntimeError("随机参数采样失败：可用点不足")

    return picked


def write_summary_txt(path, lines):
    with open(path, "w", encoding="utf-8") as f:
        for line in lines:
            f.write(line.rstrip() + "\n")


# =========================================================
# 4. 节点对比
# =========================================================

def run_node_compare():
    out_root = os.path.join(OUT_ROOT, "node_compare")
    ensure_dir(out_root)

    success = 0
    failed = 0
    summary_lines = []

    for te, tau in NODE_PARAM_SETS:
        print("=" * 72)
        print(f"[NODE] te={te:.3f}, tau={tau:.3f}")

        try:
            calspec_dir = get_calspec_dir(te, tau)
            en_file = os.path.join(calspec_dir, "en.dat")
            flux_file = os.path.join(calspec_dir, "flux_smoothed_tv.dat")

            if not os.path.exists(en_file):
                raise FileNotFoundError(f"缺少文件: {en_file}")
            if not os.path.exists(flux_file):
                raise FileNotFoundError(f"缺少文件: {flux_file}")

            E_raw = read_binary_double(en_file)
            F_raw = read_binary_double(flux_file)

            ok = np.isfinite(E_raw) & np.isfinite(F_raw) & (E_raw > 0) & (F_raw > 0)
            E_raw = E_raw[ok]
            F_raw = F_raw[ok]

            idx = np.argsort(E_raw)
            E_raw = E_raw[idx]
            F_raw = F_raw[idx]

            out_dir = os.path.join(out_root, f"te_{te:.3f}_tau_{tau:.3f}")
            ensure_dir(out_dir)

            local_qdp = os.path.join(out_dir, "local_model.qdp")
            mod_qdp   = os.path.join(out_dir, "mod_model.qdp")

            local_xcm = os.path.join(out_dir, "local_model_export.xcm") if SAVE_XCM else None
            mod_xcm   = os.path.join(out_dir, "mod_model_export.xcm") if SAVE_XCM else None

            export_local_model_qdp(te, tau, local_qdp, local_xcm)
            export_mod_model_qdp(te, tau, mod_qdp, mod_xcm)

            E_local, F_local = read_qdp_model(local_qdp)
            E_mod,   F_mod   = read_qdp_model(mod_qdp)

            F_local_i = interp_loglog(E_raw, E_local, F_local)
            F_mod_i   = interp_loglog(E_raw, E_mod, F_mod)

            F_raw_n   = normalize_at(E_raw, F_raw, E_REF)
            F_local_n = normalize_at(E_raw, F_local_i, E_REF)
            F_mod_n   = normalize_at(E_raw, F_mod_i, E_REF)

            ratio_local = F_local_n / F_raw_n
            ratio_mod   = F_mod_n / F_raw_n

            if PRINT_DEBUG:
                print_debug_block("RAW", E_raw, F_raw)
                print_debug_block("LOCAL QDP", E_local, F_local)
                print_debug_block("MOD QDP", E_mod, F_mod)
                print(f"  max |Local/TV - 1| = {np.nanmax(np.abs(ratio_local - 1.0)):.6g}")
                print(f"  max |Mod/TV   - 1| = {np.nanmax(np.abs(ratio_mod   - 1.0)):.6g}")

            out_png = os.path.join(out_dir, "node_compare.png")
            out_pdf = os.path.join(out_dir, "node_compare.pdf")

            make_node_plot(
                E_raw, F_raw_n, F_local_n, F_mod_n,
                ratio_local, ratio_mod,
                te, tau, out_png, out_pdf
            )

            summary_lines.append(
                f"NODE te={te:.3f} tau={tau:.3f} "
                f"max_abs(Local/TV-1)={np.nanmax(np.abs(ratio_local - 1.0)):.6g} "
                f"max_abs(Mod/TV-1)={np.nanmax(np.abs(ratio_mod - 1.0)):.6g}"
            )

            success += 1

        except Exception as e:
            failed += 1
            summary_lines.append(f"NODE te={te:.3f} tau={tau:.3f} FAILED: {e}")
            print(f"  FAILED: {e}")

    write_summary_txt(os.path.join(out_root, "summary.txt"), summary_lines)

    print("=" * 72)
    print(f"[NODE] 完成: 成功 {success}, 失败 {failed}")


# =========================================================
# 5. 随机插值点对比
# =========================================================

def run_random_compare():
    out_root = os.path.join(OUT_ROOT, "random_compare")
    ensure_dir(out_root)

    existing_nodes = discover_existing_nodes(BASE_EXPORT_DIR)
    random_params = choose_random_param_sets(
        existing_nodes=existing_nodes,
        n=RANDOM_COMPARE_N,
        seed=RANDOM_SEED,
        avoid_gridlike=AVOID_GRIDLIKE_POINTS
    )

    summary_lines = []
    summary_lines.append(f"RANDOM_SEED = {RANDOM_SEED}")
    summary_lines.append(f"RANDOM_COMPARE_N = {RANDOM_COMPARE_N}")
    summary_lines.append(f"TE range = [{TE_MIN}, {TE_MAX}]")
    summary_lines.append(f"TAU range = [{TAU_MIN}, {TAU_MAX}]")
    summary_lines.append("")

    success = 0
    failed = 0

    for i, (te, tau) in enumerate(random_params, start=1):
        print("=" * 72)
        print(f"[RAND {i:02d}] te={te:.3f}, tau={tau:.3f}")

        try:
            out_dir = os.path.join(out_root, f"rand_{i:02d}_te_{te:.3f}_tau_{tau:.3f}")
            ensure_dir(out_dir)

            local_qdp = os.path.join(out_dir, "local_model.qdp")
            mod_qdp   = os.path.join(out_dir, "mod_model.qdp")

            local_xcm = os.path.join(out_dir, "local_model_export.xcm") if SAVE_XCM else None
            mod_xcm   = os.path.join(out_dir, "mod_model_export.xcm") if SAVE_XCM else None

            export_local_model_qdp(te, tau, local_qdp, local_xcm)
            export_mod_model_qdp(te, tau, mod_qdp, mod_xcm)

            E_local, F_local = read_qdp_model(local_qdp)
            E_mod,   F_mod   = read_qdp_model(mod_qdp)

            # 用 local 的能量网格作为公共比较网格
            E_cmp = E_local.copy()
            F_mod_i = interp_loglog(E_cmp, E_mod, F_mod)

            F_local_n = normalize_at(E_cmp, F_local, E_REF)
            F_mod_n   = normalize_at(E_cmp, F_mod_i, E_REF)

            ratio_mod_over_local = F_mod_n / F_local_n

            if PRINT_DEBUG:
                print_debug_block("LOCAL QDP", E_local, F_local)
                print_debug_block("MOD QDP", E_mod, F_mod)
                print(f"  max |Mod/Local - 1| = {np.nanmax(np.abs(ratio_mod_over_local - 1.0)):.6g}")
                print(f"  med |Mod/Local - 1| = {np.nanmedian(np.abs(ratio_mod_over_local - 1.0)):.6g}")

            out_png = os.path.join(out_dir, "random_interp_compare.png")
            out_pdf = os.path.join(out_dir, "random_interp_compare.pdf")

            make_random_plot(
                E_cmp, F_local_n, F_mod_n, ratio_mod_over_local,
                te, tau, out_png, out_pdf
            )

            summary_lines.append(
                f"RAND {i:02d} te={te:.3f} tau={tau:.3f} "
                f"max_abs(Mod/Local-1)={np.nanmax(np.abs(ratio_mod_over_local - 1.0)):.6g} "
                f"med_abs(Mod/Local-1)={np.nanmedian(np.abs(ratio_mod_over_local - 1.0)):.6g}"
            )

            success += 1

        except Exception as e:
            failed += 1
            summary_lines.append(f"RAND {i:02d} te={te:.3f} tau={tau:.3f} FAILED: {e}")
            print(f"  FAILED: {e}")

    write_summary_txt(os.path.join(out_root, "summary.txt"), summary_lines)

    print("=" * 72)
    print(f"[RANDOM] 完成: 成功 {success}, 失败 {failed}")
    print("随机参数如下：")
    for te, tau in random_params:
        print(f"  te={te:.3f}, tau={tau:.3f}")


# =========================================================
# 6. main
# =========================================================

def main():
    if not os.path.exists(BASE_EXPORT_DIR):
        raise FileNotFoundError(f"BASE_EXPORT_DIR 不存在: {BASE_EXPORT_DIR}")
    if not os.path.exists(LOCAL_MODEL_DIR):
        raise FileNotFoundError(f"LOCAL_MODEL_DIR 不存在: {LOCAL_MODEL_DIR}")
    if not os.path.exists(MOD_MODEL_PATH):
        raise FileNotFoundError(f"MOD_MODEL_PATH 不存在: {MOD_MODEL_PATH}")

    print("=" * 72)
    print("开始节点对比")
    run_node_compare()

    print("=" * 72)
    print("开始随机插值点对比")
    run_random_compare()

    print("=" * 72)
    print("全部完成")
    print(f"输出根目录: {OUT_ROOT}")


if __name__ == "__main__":
    main()