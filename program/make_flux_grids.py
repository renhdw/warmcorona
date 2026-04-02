
"""
make_flux_grids_flat.py
- 读取“首行为能量、首列为 No”的CSV（数值为光子数谱亮度 L(E)）
- 计算光子数流量：F(E;θ,d) = L(E) * cos(θ) / (2π d^2)
- 距离 d 以 cm 计：d_cm = 10^{log10(d/Mpc)} * MPC_TO_CM
- 角度 θ：0..87 度，步长 3 度
- 距离网格：log10(d/Mpc) = 1..4，步长 0.1
- 输出为“扁平”目录：{OUTPUT_ROOT}/Dlog_<x.x>_theta_<YY>deg/spec.csv
"""

import os, json, math, csv
from typing import Tuple, List
import numpy as np
import pandas as pd

# ======== 配置（改这里） ========
INPUT_CSV   = r"/home/hdw/data/naoc/EOTA/2025.09.08/disk_spec_en.csv/disk_spec_E0.00101278-91.9984_mode-alternate_no-range-3402-10205.csv"   # 你的输入CSV
OUTPUT_ROOT = r"/home/hdw/data/naoc/EOTA/2025.09.08/disk_spec_en.csv"       # 根输出目录（会自动创建）

ANGLE_START, ANGLE_END, ANGLE_STEP = 0, 87, 30     # 角度设置（度）
LOGD_START, LOGD_END, LOGD_STEP    = 1.0, 4.0, 2# log10(d/Mpc) 设置

COERCE_NUMERIC = False  # 如需将有效能量列强制为数值，可置 True
# =================================

MPC_TO_CM = 3.08567758e24

def read_csv_loose(path: str) -> pd.DataFrame:
    return pd.read_csv(path, header=None, sep=None, engine="python")

def parse_energy_row(energy_row: pd.Series) -> Tuple[List[int], np.ndarray]:
    valid_cols, energies = [], []
    for col in range(1, len(energy_row)):  # 跳过 No 列
        try:
            e = float(energy_row.iloc[col])
            if math.isfinite(e):
                valid_cols.append(col)
                energies.append(e)
        except Exception:
            pass
    if not valid_cols:
        raise ValueError("首行未解析到任何有效能量列。请确认首行是能量数字。")
    return valid_cols, np.array(energies, dtype=float)

def build_angle_list(start: int, end: int, step: int) -> np.ndarray:
    return np.arange(start, end + 1e-9, step, dtype=int)

def build_logd_list(start: float, end: float, step: float) -> np.ndarray:
    n = int(round((end - start) / step)) + 1
    vals = start + np.arange(n, dtype=float) * step
    return np.round(vals, 1)

def scale_df_to_flux(df_in: pd.DataFrame,
                     valid_cols: List[int],
                     scale: float,
                     coerce_numeric: bool = False) -> pd.DataFrame:
    df = df_in.copy(deep=True)
    if coerce_numeric:
        for c in valid_cols:
            df.iloc[1:, c] = pd.to_numeric(df.iloc[1:, c], errors="coerce")
    df.iloc[1:, valid_cols] = df.iloc[1:, valid_cols] * scale
    return df

def write_spec_csv(df_out: pd.DataFrame, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "spec.csv")
    df_out.to_csv(out_csv, header=False, index=False, na_rep="", quoting=csv.QUOTE_MINIMAL)
    return out_csv

def main():
    df = read_csv_loose(INPUT_CSV)
    valid_cols, energies = parse_energy_row(df.iloc[0, :])

    angle_list = build_angle_list(ANGLE_START, ANGLE_END, ANGLE_STEP)   # 度
    logd_list  = build_logd_list(LOGD_START, LOGD_END, LOGD_STEP)       # 无量纲 log10(d/Mpc)

    two_pi = 2.0 * math.pi
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    written = []
    for logd in logd_list:
        d_cm = (10.0 ** logd) * MPC_TO_CM
        for theta_deg in angle_list:
            theta_rad = math.radians(float(theta_deg))
            scale = math.cos(theta_rad) / (two_pi * (d_cm ** 2))

            df_flux = scale_df_to_flux(df, valid_cols, scale, coerce_numeric=COERCE_NUMERIC)

            # 扁平目录名：Dlog_<x.x>_theta_<YY>deg
            out_dir = os.path.join(OUTPUT_ROOT, f"Dlog_{logd:.1f}_theta_{theta_deg:02d}deg")
            out_csv = write_spec_csv(df_flux, out_dir)
            written.append({
                "out_csv": out_csv,
                "logd": float(logd),
                "theta_deg": int(theta_deg),
                "d_cm": d_cm,
                "scale": scale
            })

    # 写运行日志
    run_log = {
        "input_csv": os.path.abspath(INPUT_CSV),
        "output_root": os.path.abspath(OUTPUT_ROOT),
        "angle_deg": {
            "start": ANGLE_START, "end": ANGLE_END, "step": ANGLE_STEP,
            "count": int(len(angle_list))
        },
        "log10_d_Mpc": {
            "start": LOGD_START, "end": LOGD_END, "step": LOGD_STEP,
            "count": int(len(logd_list))
        },
        "unit_notes": {
            "distance": "d is converted to cm: d_cm = 10^{log10(d/Mpc)} * 3.08567758e24 cm",
            "formula": "F(E;theta,d) = L(E) * cos(theta) / (2*pi*d^2)",
            "flux_unit_hint": "If L(E) is photons/s/keV, F is photons/cm^2/s/keV"
        },
        "energies_keV_count": int(len(energies)),
        "example_outputs": [w["out_csv"] for w in written[:5]],
        "total_files": len(written)
    }
    with open(os.path.join(OUTPUT_ROOT, "run.log"), "w", encoding="utf-8") as f:
        f.write(json.dumps(run_log, ensure_ascii=False, indent=2))

    print(f"✅ 完成，共写出 {len(written)} 个 spec.csv")
    if written:
        print("示例：", written[0]["out_csv"])
    print(f"运行日志：{os.path.join(OUTPUT_ROOT, 'run.log')}")

if __name__ == "__main__":
    main()
