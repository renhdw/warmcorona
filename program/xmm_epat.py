#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模块 X：PN 模式诊断（epatplot 专用）

前提：
  观测目录下已经有：
    - PN.fits         # 模块2 epproc 生成的 PN 事件
    - PNgti.fits      # 模块3 tabgtigen 生成的 GTI
    - ccf.cif         # 模块2 cifbuild 生成
    - *SUM.SAS        # 模块2 odfingest 生成
    - src.reg         # DS9 画的源区（使用 X,Y 坐标）

目录结构：
  ~/data/XMM/<SOURCE>/<OBSID>/
      PN.fits
      PNgti.fits
      ccf.cif
      ...SUM.SAS
      src.reg
      （本脚本生成）
      pn_filtered.evt
      pn_filtered_pat.ps

功能：
  1. 设置 SAS_CCF / SAS_ODF / ODF_PATH
  2. 从 src.reg 解析 circle(x,y,r)
  3. evselect 过滤：区域 + GTI → pn_filtered.evt
  4. 设置 MPLBACKEND=Agg（避免调 X11）
  5. epatplot set=pn_filtered.evt plotfile=pn_filtered_pat.ps
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd, cwd: Path):
    """在 cwd 下运行命令，并把输出直接打印出来。"""
    print(f"[CMD@{cwd.name}] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(cwd))
    if result.returncode != 0:
        raise RuntimeError(f"命令失败：{' '.join(cmd)} (returncode={result.returncode})")


def setup_sas_env(obs_dir: Path):
    """
    设置 SAS 环境变量：
      SAS_CCF = ./ccf.cif
      SAS_ODF = ./<something>SUM.SAS
      ODF_PATH = obs_dir
    """
    ccf = obs_dir / "ccf.cif"
    if not ccf.exists():
        raise RuntimeError(f"[env] 未找到 ccf.cif，请确认在 {obs_dir} 跑过 cifbuild。")

    sas_files = list(obs_dir.glob("*SUM.SAS"))
    if not sas_files:
        raise RuntimeError(f"[env] 未在 {obs_dir} 找到 *SUM.SAS，请确认跑过 odfingest。")
    if len(sas_files) > 1:
        print(f"[env] 警告：找到多个 SUM.SAS，将使用：{sas_files[0].name}")
    sum_sas = sas_files[0]

    os.environ["SAS_CCF"] = str(ccf)
    os.environ["SAS_ODF"] = str(sum_sas)
    os.environ["ODF_PATH"] = str(obs_dir)

    print("[env] SAS_CCF =", os.environ["SAS_CCF"])
    print("[env] SAS_ODF =", os.environ["SAS_ODF"])
    print("[env] ODF_PATH =", os.environ["ODF_PATH"])


def parse_ds9_region(reg_file: Path):
    """
    解析 DS9 region 文件中的 circle(x,y,r)，返回 (x, y, r).

    只支持单个圆：
        circle(27173.154,27149.598,990.544)
    """
    if not reg_file.exists():
        raise RuntimeError(f"region 文件不存在：{reg_file}")
    txt = reg_file.read_text()
    m = re.search(r"circle\(([^,]+),([^,]+),([^)]+)\)", txt)
    if not m:
        raise RuntimeError(f"无法在 {reg_file} 中解析 circle(x,y,r)，请确认是 DS9 circle 区域。")
    x = float(m.group(1))
    y = float(m.group(2))
    r = float(m.group(3))
    return x, y, r


def run_epat_for_obsid(obs_dir: Path, reg_name: str = "src.reg",
                       pn_name: str = "PN.fits", gti_name: str = "PNgti.fits"):
    """
    对单个 OBSID 目录执行：
      - 设置 SAS 环境变量
      - 从 src.reg 解析 circle(x,y,r)
      - evselect → pn_filtered.evt
      - 设置 MPLBACKEND=Agg
      - epatplot → pn_filtered_pat.ps
    """
    print("=" * 60)
    print(f"[epatplot 模块] 观测目录：{obs_dir}")
    print("=" * 60)

    setup_sas_env(obs_dir)

    pn_evt = obs_dir / pn_name
    gti = obs_dir / gti_name
    reg_file = obs_dir / reg_name

    if not pn_evt.exists():
        raise RuntimeError(f"未找到 PN 事件文件：{pn_evt}")
    if not gti.exists():
        raise RuntimeError(f"未找到 GTI 文件：{gti}")
    if not reg_file.exists():
        raise RuntimeError(f"未找到 region 文件：{reg_file}")

    print(f"[info] 使用事件文件：{pn_evt.name}")
    print(f"[info] 使用 GTI：{gti.name}")
    print(f"[info] 使用区域文件：{reg_file.name}")

    # 1) 解析 region
    x, y, r = parse_ds9_region(reg_file)
    print(f"[region] circle({x:.3f}, {y:.3f}, {r:.3f})")

    # 2) evselect 过滤到 pn_filtered.evt
    filtered = obs_dir / "pn_filtered.evt"
    expr = f"((X,Y) IN CIRCLE({x},{y},{r})) && gti({gti.name},TIME)"
    print(f"[evselect] expression = {expr}")
    cmd = [
        "evselect",
        f"table={pn_evt.name}",
        "withfilteredset=yes",
        f"filteredset={filtered.name}",
        "keepfilteroutput=yes",
        f"expression={expr}",
    ]
    run_cmd(cmd, obs_dir)
    if not filtered.exists():
        raise RuntimeError("evselect 未生成 pn_filtered.evt")

    print(f"[evselect] 已生成过滤事件文件：{filtered}")

    # 3) 设置 MPLBACKEND=Agg（避免调 X11/matplotlib GUI）
    os.environ["MPLBACKEND"] = "Agg"
    print("[env] MPLBACKEND=Agg 已设置。")

    # 4) epatplot
    pat_ps = obs_dir / "pn_filtered_pat.ps"
    cmd = [
        "epatplot",
        f"set={filtered.name}",
        f"plotfile={pat_ps.name}",
    ]
    run_cmd(cmd, obs_dir)

    print(f"[epatplot] 已生成模式诊断图：{pat_ps}")
    print("=" * 60)
    print()


def main():
    parser = argparse.ArgumentParser(
        description="PN 模式诊断：PN.fits + PNgti.fits + src.reg → pn_filtered.evt + epatplot"
    )
    parser.add_argument("source", help="源名，例如 1H0419-577")
    parser.add_argument(
        "obsids",
        nargs="+",
        help="观测号（一个或多个），例如：0148000201 0148000301",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "data" / "XMM",
        help="数据根路径（默认 ~/data/XMM）",
    )
    parser.add_argument(
        "--reg-name",
        default="src.reg",
        help="region 文件名（默认 src.reg）",
    )
    parser.add_argument(
        "--pn-name",
        default="PN.fits",
        help="PN 事件文件名（默认 PN.fits）",
    )
    parser.add_argument(
        "--gti-name",
        default="PNgti.fits",
        help="GTI 文件名（默认 PNgti.fits）",
    )

    args = parser.parse_args()
    data_root: Path = args.data_root.expanduser()

    for obsid in args.obsids:
        obs_dir = data_root / args.source / obsid
        if not obs_dir.exists():
            print(f"[错误] 目录不存在，跳过：{obs_dir}", file=sys.stderr)
            continue

        try:
            run_epat_for_obsid(
                obs_dir,
                reg_name=args.reg_name,
                pn_name=args.pn_name,
                gti_name=args.gti_name,
            )
        except Exception as e:
            print(f"[错误] 处理 {obsid} 时失败：{e}", file=sys.stderr)


if __name__ == "__main__":
    main()
