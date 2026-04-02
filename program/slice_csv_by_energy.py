#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
slice_csv_by_energy.py (简洁版)
- 不做任何数值转换或dtype检查
- 功能: 能量列范围截取 + 列模式(all/alternate) + No行筛选(all/range/list)
- 输出命名基于“实际选中能量范围 + 列模式 + No规则”，并写log
"""

import os, math, json, csv, re
from typing import List, Tuple
import numpy as np
import pandas as pd

# ====== 你可直接修改 ======
INPUT_CSV  = r"/home/hdw/data/naoc/EOTA/2025.09.08/disk_spec.csv"
OUTPUT_DIR = r"/home/hdw/data/naoc/EOTA/2025.09.08/disk_spec_en.csv"  # 当作目录，会自动创建（即使带.csv后缀也可作为目录名）

EMIN, EMAX = 0.001, 100.0          # keV 请求范围
DEFAULT_MODE = "all"               # "all" 或 "alternate"

NO_FILTER_MODE = "all"             # "all" / "range" / "list"
NO_RANGES = [(0, 10)]              # range时使用(闭区间，可多段)
NO_LIST   = [1, 2, 3]              # list时使用

INTERACTIVE_PICK_MODE = True       # 运行时可手选 all/alternate
INTERACTIVE_NO_FILTER = True       # 运行时可手选 No 的规则
# ==========================

def read_csv_loose(path: str) -> pd.DataFrame:
    """首行是能量，首列是No，不设表头。"""
    return pd.read_csv(path, header=None, sep=None, engine="python")

def parse_energy_row(energy_row: pd.Series) -> Tuple[List[int], np.ndarray]:
    """解析首行能量（跳过第0列No），返回有效列索引和能量数组。"""
    valid_cols, energies = [], []
    for col in range(1, len(energy_row)):
        try:
            e = float(energy_row.iloc[col])
            if math.isfinite(e):
                valid_cols.append(col)
                energies.append(e)
        except Exception:
            pass
    if not valid_cols:
        raise ValueError("首行未解析到任何有效能量列。")
    return valid_cols, np.array(energies, dtype=float)

def nearest_range(energies: np.ndarray, emin: float, emax: float) -> Tuple[int, int]:
    """找到最接近 [emin, emax] 的两个能量索引（按绝对差最小），返回升序元组。"""
    i_lo = int(np.argmin(np.abs(energies - emin)))
    i_hi = int(np.argmin(np.abs(energies - emax)))
    return tuple(sorted([i_lo, i_hi]))

def pick_columns(valid_df_cols: List[int], lo: int, hi: int, mode: str) -> List[int]:
    """按模式选择列；总是保留第0列No。"""
    cols = valid_df_cols[lo:hi+1]
    if mode == "alternate":
        cols = cols[::2]
    return [0] + cols

def build_no_mask(df: pd.DataFrame, mode: str,
                  ranges: List[Tuple[float, float]],
                  id_list: List[float]) -> np.ndarray:
    """根据 No 选择行；首行(能量行)始终保留。"""
    n = len(df)
    mask = np.zeros(n, dtype=bool)
    mask[0] = True  # 能量行
    if mode == "all":
        mask[1:] = True
        return mask

    no_vals = pd.to_numeric(df.iloc[1:, 0], errors="coerce").to_numpy()
    ok = np.isfinite(no_vals)

    if mode == "range":
        rowmask = np.zeros(len(no_vals), dtype=bool)
        for lo, hi in ranges:
            if lo > hi: lo, hi = hi, lo
            rowmask |= (no_vals >= lo) & (no_vals <= hi)
        mask[1:] = ok & rowmask
    elif mode == "list":
        s = set(id_list)
        rowmask = np.array([v in s for v in no_vals])
        mask[1:] = ok & rowmask
    else:
        raise ValueError("NO_FILTER_MODE 必须是 all/range/list")
    return mask

def safe_tag(s: str) -> str:
    """文件名安全化。"""
    return re.sub(r"[^A-Za-z0-9._-]+", "-", s).strip("-")

def build_no_tag(mode: str, ranges: List[Tuple[float, float]], ids: List[float]) -> str:
    """构建 No 选择的标签文本，用于文件名。"""
    if mode == "all":
        return "no-all"
    if mode == "range":
        segs = []
        for lo, hi in ranges:
            if lo > hi: lo, hi = hi, lo
            segs.append(f"{lo:g}-{hi:g}")
        return "no-range-" + "_".join(safe_tag(x) for x in segs) if segs else "no-range-empty"
    if mode == "list":
        segs = [f"{v:g}" for v in ids]
        return "no-list-" + "_".join(safe_tag(x) for x in segs) if segs else "no-list-empty"
    return "no-unknown"

def save_with_log(df_in: pd.DataFrame,
                  row_mask: np.ndarray,
                  keep_cols: List[int],
                  energies_all: np.ndarray,
                  valid_cols: List[int],
                  emin_req: float, emax_req: float,
                  mode: str,
                  input_path: str,
                  output_dir: str,
                  no_filter_conf: dict) -> Tuple[str, str]:
    """
    保存筛选后的 CSV 和日志。注意：这里的 no_filter_conf 是字典，包含
      {'mode': ..., 'ranges': ..., 'list': ...}
    """
    os.makedirs(output_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(input_path))[0]

    # 行、列筛选（不改dtype/精度）
    out_df = df_in.loc[row_mask, keep_cols]

    # 实际能量范围（按选中列计算）
    kept_energy_cols = [c for c in keep_cols if c in valid_cols]
    kept_idx = [valid_cols.index(c) for c in kept_energy_cols]
    energies_kept = [float(energies_all[i]) for i in kept_idx]
    if energies_kept:
        e_min, e_max = min(energies_kept), max(energies_kept)
        e_tag = f"E{e_min:g}-{e_max:g}"
    else:
        e_tag = "E-none"

    base = safe_tag(f"{stem}_{e_tag}_mode-{mode}_"
                    f"{build_no_tag(no_filter_conf['mode'], no_filter_conf['ranges'], no_filter_conf['list'])}")
    out_csv = os.path.join(output_dir, base + ".csv")
    out_log = os.path.join(output_dir, base + ".log")

    # 写CSV：不指定 float_format，保持 pandas 默认（不强行截断位数）
    out_df.to_csv(out_csv, header=False, index=False, na_rep="", quoting=csv.QUOTE_MINIMAL)

    # 记录 No（去掉能量行）
    kept_no = pd.to_numeric(out_df.iloc[1:, 0], errors="coerce").dropna().tolist() if out_df.shape[0] > 1 else []

    # 写log
    log = {
        "input_csv": os.path.abspath(input_path),
        "output_csv": os.path.abspath(out_csv),
        "energy_request_keV": {"emin": emin_req, "emax": emax_req},
        "energy_kept_keV": {
            "min": (min(energies_kept) if energies_kept else None),
            "max": (max(energies_kept) if energies_kept else None),
            "count": len(energies_kept)
        },
        "column_pick_mode": mode,
        "kept_dataframe_columns": keep_cols,
        "no_filter": no_filter_conf,
        "kept_no_values": kept_no
    }
    with open(out_log, "w", encoding="utf-8") as f:
        f.write(json.dumps(log, ensure_ascii=False, indent=2))

    return out_csv, out_log

def main():
    df = read_csv_loose(INPUT_CSV)

    # 解析能量、匹配范围
    valid_cols, energies = parse_energy_row(df.iloc[0, :])
    lo, hi = nearest_range(energies, EMIN, EMAX)

    # 列模式（可交互）
    mode = DEFAULT_MODE
    if INTERACTIVE_PICK_MODE:
        sel = input("选列模式 [all/alternate]（回车默认 all）: ").strip().lower()
        if sel in ("all", "alternate"):
            mode = sel
        elif sel != "":
            print("输入无效，使用默认 all")
            mode = "all"
    keep_cols = pick_columns(valid_cols, lo, hi, mode)

    # No 行筛选（可交互）
    nf_mode, ranges, id_list = NO_FILTER_MODE, NO_RANGES, NO_LIST
    if INTERACTIVE_NO_FILTER:
        sel = input("No 筛选方式 [all/range/list]（回车默认 all）: ").strip().lower()
        if sel in ("all", "range", "list"):
            nf_mode = sel
        elif sel != "":
            print("输入无效，使用默认 all")
            nf_mode = "all"

        if nf_mode == "range":
            txt = input("输入区间，如 0-100, 500-800（闭区间）: ").strip()
            parsed = []
            if txt:
                for seg in txt.split(","):
                    if "-" in seg:
                        a, b = seg.split("-", 1)
                        try:
                            lo2, hi2 = float(a.strip()), float(b.strip())
                            if lo2 > hi2:
                                lo2, hi2 = hi2, lo2
                            parsed.append((lo2, hi2))
                        except Exception:
                            pass
            if parsed:
                ranges = parsed
            else:
                print("未解析到有效区间，退回全选。")
                nf_mode = "all"

        elif nf_mode == "list":
            txt = input("输入编号列表，如 1,2,3,10: ").strip()
            parsed = []
            if txt:
                for x in txt.split(","):
                    try:
                        parsed.append(float(x.strip()))
                    except Exception:
                        pass
            if parsed:
                id_list = parsed
            else:
                print("未解析到有效编号，退回全选。")
                nf_mode = "all"

    row_mask = build_no_mask(df, nf_mode, ranges, id_list)
    no_filter_conf = {"mode": nf_mode, "ranges": ranges, "list": id_list}

    out_csv, out_log = save_with_log(
        df_in=df, row_mask=row_mask, keep_cols=keep_cols,
        energies_all=energies, valid_cols=valid_cols,
        emin_req=EMIN, emax_req=EMAX, mode=mode,
        input_path=INPUT_CSV, output_dir=OUTPUT_DIR,
        no_filter_conf=no_filter_conf
    )

    print("✅ 完成")
    print("CSV ->", out_csv)
    print("LOG ->", out_log)

if __name__ == "__main__":
    main()
