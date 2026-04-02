#!/usr/bin/env python3
"""
模块 2 + 模块 3 ：XMM-SAS 前期处理 + 耀斑背景筛选（仅 PN）

模块 2：前期数据处理
    对每个 OBSID 目录执行：
      1. 设置环境变量 (SAS_ODF / SAS_CCF)
      2. cifbuild → ccf.cif
      3. odfingest → *SUM.SAS
      4. 设置 SAS_ODF 为 SUM.SAS
      5. epproc → 生成 PN 事件文件
      6. 将 PN 事件文件重命名为 PN.fits

模块 3：耀斑背景处理（flare filtering）
    在每个 OBSID 目录中，对 PN.fits 执行：
      1. 用高能段 (10–12 keV, PATTERN==0) 抽 light curve → ratePN.fits
      2. 手动可用 dsplot 检查 flaring (脚本只打印提示)
      3. 使用 tabgtigen 生成 GTI：PNgti.fits
      4. 用 GTI 清洗 PN.fits → PNclean.fits

运行示例：
    python3 xmm_preprocess.py 1H0419-577 0148000301
    python3 xmm_preprocess.py 1H0419-577 0148000201 0148000301
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


# ============================================
# 通用工具函数
# ============================================

def run_cmd(cmd: list[str], cwd: Path):
    """
    在给定工作目录 cwd 下执行命令 cmd（列表形式），
    并将 stdout/stderr 直接打印到当前终端。
    """
    print(f"[CMD@{cwd.name}] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(cwd))
    if result.returncode != 0:
        raise RuntimeError(f"命令失败：{' '.join(cmd)} (returncode={result.returncode})")


def find_sum_sas(obs_dir: Path) -> Path:
    """
    在观测目录内查找唯一的 *SUM.SAS 文件。
    这是 odfingest 的输出，也是后续 SAS_ODF 应指向的文件。
    """
    sas_files = list(obs_dir.glob("*SUM.SAS"))
    if not sas_files:
        raise RuntimeError(f"未找到 SUM.SAS 文件：{obs_dir}")
    if len(sas_files) > 1:
        print(f"[警告] 找到多个 SUM.SAS，将使用第一个：{sas_files[0].name}")
    return sas_files[0]


def find_pn_event(obs_dir: Path) -> Path:
    """
    查找 epproc 生成的 PN 事件文件（原始名字比较长）。
    通常形如：0045_0120300201_EPN_U002_ImagingEvts.ds / .FIT
    """
    candidates = list(obs_dir.glob("*EPN*Evts*.ds")) + list(obs_dir.glob("*EPN*Evts*.FIT"))
    if not candidates:
        raise RuntimeError(f"epproc 未生成 PN 事件文件：{obs_dir}")
    return candidates[0]


# ============================================
# 模块 2：SAS 前期处理
# ============================================

def run_sas_initial(obs_dir: Path) -> Path:
    """
    模块 2：在单个 OBSID 目录中执行：
      - 设置环境变量
      - cifbuild
      - odfingest
      - 设置 SAS_ODF
      - epproc
      - 将 PN 文件重命名为 PN.fits

    返回值：
      PN 事件文件路径：obs_dir / "PN.fits"
    """
    print("=" * 60)
    print(f"[模块2] 前期处理：{obs_dir}")
    print("=" * 60)

    cwd = obs_dir

    # 1. 设置环境变量
    print("[2-1] 设置 SAS 环境变量")
    os.environ["ODF_PATH"] = str(obs_dir)   # 方便你自己记
    os.environ["SAS_ODF"] = str(obs_dir)    # 初始先指向目录
    os.environ["SAS_CCF"] = str(obs_dir)    # cifbuild 会在此产生 ccf.cif

    # 2. cifbuild → ccf.cif
    print("[2-2] 运行 cifbuild → 生成 ccf.cif")
    run_cmd(["cifbuild"], cwd)
    cif = cwd / "ccf.cif"
    if not cif.exists():
        raise RuntimeError("cifbuild 未产生 ccf.cif")
    os.environ["SAS_CCF"] = str(cif)
    print(f"    SAS_CCF = {cif}")

    # 3. odfingest → *SUM.SAS
    print("[2-3] 运行 odfingest → 生成 SUM.SAS")
    run_cmd(["odfingest"], cwd)
    sum_sas = find_sum_sas(cwd)

    # 4. 设置 SAS_ODF 指向 SUM.SAS
    os.environ["SAS_ODF"] = str(sum_sas)
    print(f"    SAS_ODF = {sum_sas}")

    # 5. epproc → PN 事件文件
    print("[2-4] 运行 epproc（PN 处理）")
    run_cmd(["epproc"], cwd)

    # 6. 找到 PN 事件文件并重命名为 PN.fits
    pn_evt = find_pn_event(cwd)
    final_pn = cwd / "PN.fits"
    print(f"[2-5] 找到 PN 事件文件：{pn_evt.name} → 重命名为 PN.fits")
    if final_pn.exists():
        print("    警告：PN.fits 已存在，将被覆盖")
        final_pn.unlink()
    pn_evt.rename(final_pn)

    print(f"[模块2完成] 生成 PN.fits：{final_pn}")
    print()
    return final_pn


# ============================================
# 模块 3：耀斑背景处理（flare filtering）
# ============================================

def pn_make_high_energy_lc(obs_dir: Path, pn_event: Path) -> Path:
    """
    从 PN 事件表生成高能 light curve：
      - 输入：PN.fits
      - 输出：ratePN.fits（高能 LC，10–12 keV, PATTERN==0）

    对应命令：
      evselect table=PN.fits withrateset=Y rateset=ratePN.fits \
               maketimecolumn=Y timebinsize=100 makeratecolumn=Y \
               expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'
    """
    print("[3-1] 生成高能 LC：ratePN.fits")
    cwd = obs_dir
    rate_file = cwd / "ratePN.fits"

    cmd = [
        "evselect",
        f"table={pn_event.name}",
        "withrateset=Y",
        f"rateset={rate_file.name}",
        "maketimecolumn=Y",
        "timebinsize=100",
        "makeratecolumn=Y",
        "expression=#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)",
    ]
    run_cmd(cmd, cwd)
    if not rate_file.exists():
        raise RuntimeError("未生成 ratePN.fits")
    print(f"    生成：{rate_file}")
    return rate_file


def pn_generate_gti_from_rate(obs_dir: Path, rate_file: Path, rate_cut: float = 0.4) -> Path:
    """
    使用 tabgtigen 从高能 LC 生成 GTI：
      - 输入：ratePN.fits
      - 输出：PNgti.fits

    对应命令：
      tabgtigen table=ratePN.fits expression='RATE<=0.4' gtiset=PNgti.fits
    """
    print("[3-2] 使用 tabgtigen 生成 GTI：PNgti.fits")
    cwd = obs_dir
    gti_file = cwd / "PNgti.fits"

    expr = f"RATE<={rate_cut}"
    cmd = [
        "tabgtigen",
        f"table={rate_file.name}",
        f"expression={expr}",
        f"gtiset={gti_file.name}",
    ]
    run_cmd(cmd, cwd)
    if not gti_file.exists():
        raise RuntimeError("未生成 GTI：PNgti.fits")
    print(f"    生成：{gti_file}")
    return gti_file


def pn_apply_gti_clean_events(obs_dir: Path, pn_event: Path, gti_file: Path) -> Path:
    """
    使用 GTI 文件清洗 PN 事件表，生成 PNclean.fits：

    对应命令：
      evselect table=PN.fits withfilteredset=Y filteredset=PNclean.fits \
               destruct=Y keepfilteroutput=T \
               expression='#XMMEA_EP && gti(PNgti.fits,TIME) && (PI>150)'
    """
    print("[3-3] 应用 GTI 清洗 PN 事件 → PNclean.fits")
    cwd = obs_dir
    clean_file = cwd / "PNclean.fits"

    cmd = [
        "evselect",
        f"table={pn_event.name}",
        "withfilteredset=Y",
        f"filteredset={clean_file.name}",
        "destruct=Y",
        "keepfilteroutput=T",
        "expression=#XMMEA_EP && gti(PNgti.fits,TIME) && (PI>150)",
    ]
    run_cmd(cmd, cwd)
    if not clean_file.exists():
        raise RuntimeError("未生成 PNclean.fits")
    print(f"    生成：{clean_file}")
    return clean_file


def run_flare_filtering(obs_dir: Path, pn_event: Path, rate_cut: float = 0.4):
    """
    模块 3 总控函数（仅 PN）：

      1. 生成高能 LC：ratePN.fits
      2. 提示你可以用 dsplot 查看：
           dsplot table=ratePN.fits x=TIME y=RATE.ERROR
      3. 用 tabgtigen 生成 PNgti.fits
      4. 用 PNgti 清洗 PN.fits → PNclean.fits
    """
    print("=" * 60)
    print(f"[模块3] 耀斑背景处理（PN）：{obs_dir}")
    print("=" * 60)

    # 3-1: 高能 LC
    rate_file = pn_make_high_energy_lc(obs_dir, pn_event)

    print("    可选手动检查（不在脚本中执行）：")
    print("      dsplot table=ratePN.fits x=TIME y=RATE.ERROR")

    # 3-2: 生成 GTI
    gti_file = pn_generate_gti_from_rate(obs_dir, rate_file, rate_cut=rate_cut)

    # 3-3: 清洗 PN 事件
    clean_file = pn_apply_gti_clean_events(obs_dir, pn_event, gti_file)

    print(f"[模块3完成] PNclean.fits 已生成：{clean_file}")
    print()
    return clean_file


# ============================================
# CLI 主入口
# ============================================

def main():
    parser = argparse.ArgumentParser(
        description="模块2+3：XMM-SAS 前期处理 + 耀斑筛选（PN）。"
    )
    parser.add_argument("source", help="源名，例如：1H0419-577")
    parser.add_argument("obsids", nargs="+", help="观测号 OBSID，例如：0148000301")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "data" / "XMM",
        help="数据根目录（默认：~/data/XMM）",
    )
    parser.add_argument(
        "--rate-cut",
        type=float,
        default=0.4,
        help="tabgtigen 的 RATE 阈值（默认 0.4）",
    )
    args = parser.parse_args()

    data_root: Path = args.data_root.expanduser()

    for obsid in args.obsids:
        obs_dir = data_root / args.source / obsid
        if not obs_dir.exists():
            print(f"[错误] 找不到观测目录：{obs_dir}（请先用模块1下载+解压）")
            continue

        try:
            # 模块 2：SAS 初始处理
            pn_event = run_sas_initial(obs_dir)

            # 模块 3：耀斑背景筛选（PN）
            run_flare_filtering(obs_dir, pn_event, rate_cut=args.rate_cut)

        except Exception as e:
            print(f"[错误] 处理 {obsid} 时失败：{e}", file=sys.stderr)


if __name__ == "__main__":
    main()
