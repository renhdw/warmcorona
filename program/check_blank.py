#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

# 文件路径
csv_path = "/home/hdw/data/naoc/EOTA/2025.09.08/disk_spec.csv"

# 读入（无表头）
df = pd.read_csv(csv_path, header=None)

# 取前四列
df4 = df.iloc[:, :4]

# 判断空白：NaN 或 空字符串
is_blank = df4.isna() | df4.applymap(lambda x: str(x).strip() == "")

# 输出结果
if is_blank.values.any():
    print("⚠️ 前四列存在空白数据：")
    rows, cols = is_blank.to_numpy().nonzero()
    for r, c in zip(rows, cols):
        print(f"第 {r+1} 行, 第 {c+1} 列为空白")
else:
    print("✅ 前四列没有空白数据")
