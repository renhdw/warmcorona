#!/usr/bin/env python3
"""
模块 4：PN 光谱提取（手动选区，一条龙，交互式可选 OOT 校正）

前提：
  目标目录下已经有：
    - PNclean.fits        # 模块2+3生成的清洗后的 PN 事件表
    - ccf.cif             # 模块2 cifbuild 生成
    - *SUM.SAS            # 模块2 odfingest 生成

目录结构约定：
  ~/data/XMM/<SOURCE>/<OBSID>/
      PNclean.fits
      ccf.cif
      ...SUM.SAS
      （运行后生成：）
      PNimage.fits
      src.reg
      bg.reg
      PNsource_spectrum.fits
      PNbackground_spectrum.fits
      （如做 OOT：）
      PNsource_OOT_spectrum.fits
      PNbackground_OOT_spectrum.fits
      （注意：PNsource_spectrum.fits / PNbackground_spectrum.fits
              会被就地改成减完 OOT 的版本）
      PN.rmf
      PN.arf
      PN_spectrum_grp.fits
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


# ================================
# 公共工具函数
# ================================

def run_cmd(cmd, cwd: Path):
    """
    在 cwd 下运行命令，并将 stdout/stderr 直接打印到当前终端。
    若命令返回非 0，抛出 RuntimeError。
    """
    print(f"[CMD@{cwd.name}] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(cwd))
    if result.returncode != 0:
        raise RuntimeError(f"命令失败：{' '.join(cmd)}")


def setup_sas_env(obs_dir: Path):
    """
    设置 SAS 环境变量：
      - SAS_CCF = ./ccf.cif
      - SAS_ODF = ./<something>SUM.SAS
      - ODF_PATH = obs_dir
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

    只支持形如：
        circle(26088.8,27992.9,40)
    的单一圆区域。
    """
    txt = reg_file.read_text()
    m = re.search(r"circle\(([^,]+),([^,]+),([^)]+)\)", txt)
    if not m:
        raise RuntimeError(f"无法在 {reg_file} 中解析 circle(x,y,r)，请确认是 DS9 circle 区域。")
    x = float(m.group(1))
    y = float(m.group(2))
    r = float(m.group(3))
    return x, y, r


# ================================
# 图像生成 + 手动画圈
# ================================

def ensure_image_and_regions(obs_dir: Path, pn_clean: Path):
    """
    1. 如无 PNimage.fits，使用 PNclean.fits 生成成像（X,Y 空间）。
    2. 如无 src.reg / bg.reg：
         - 自动启动 ds9 <PNimage.fits> &
         - 提示用户在 ds9 中画 src.reg / bg.reg
    3. 检查 src.reg / bg.reg 存在后返回它们的路径。
    """
    cwd = obs_dir
    image_file = obs_dir / "PNimage.fits"

    # 1) 若无图像，则生成
    if not image_file.exists():
        print("[4-1] 未找到 PNimage.fits，开始生成成像用于选区...")
        cmd = [
            "evselect",
            f"table={pn_clean.name}",
            "imagebinning=binSize",
            f"imageset={image_file.name}",
            "withimageset=yes",
            "xcolumn=X", "ycolumn=Y",
            "ximagebinsize=80",
            "yimagebinsize=80",
        ]
        run_cmd(cmd, cwd)
        print(f"      已生成图像：{image_file}")
    else:
        print(f"[4-1] 已存在 PNimage.fits：{image_file}，跳过生成。")

    src_reg = obs_dir / "src.reg"
    bg_reg  = obs_dir / "bg.reg"

    # 2) 若 region 文件不存在，则自动启动 ds9，并提示手动画圈
    if not (src_reg.exists() and bg_reg.exists()):
        print("\n🟦 请选择源区和背景区（手动模式）")
        print("  当前工作目录：", obs_dir)
        print("  自动尝试启动 DS9 显示 PNimage.fits ...")

        try:
            subprocess.Popen(
                ["ds9", str(image_file)],
                cwd=str(cwd),
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            print(f"  ✅ 已调用 ds9 {image_file}")
        except FileNotFoundError:
            print("  ⚠️ 未找到 ds9 可执行文件，请确认已安装并在 PATH 中。")
            print(f"     你可以手动运行：cd {obs_dir} && ds9 PNimage.fits &")

        print("\n  操作说明：")
        print("    1) 在 DS9 中，用 Circle 工具画出源区域：")
        print("         - 建议：先画略大的圈包住源 → Region → Centroid → 调整半径")
        print(f"       然后在 Region 窗口中保存为：{src_reg.name}")
        print("         （目录就是当前 obs 目录）")
        print("    2) 再画背景区域（同一 CCD，远离源，无点源，半径与源相近），")
        print(f"       保存为：{bg_reg.name}")
        input("\n画好两个圈并保存之后，回到终端按回车继续...")

    # 3) 再次检查 region 文件是否存在
    if not src_reg.exists():
        raise RuntimeError(f"未找到 src.reg：{src_reg}")
    if not bg_reg.exists():
        raise RuntimeError(f"未找到 bg.reg：{bg_reg}")

    print(f"[4-2] 已检测到 src.reg 和 bg.reg，可继续光谱提取。")
    return src_reg, bg_reg


# ================================
# OOT 相关函数（按你笔记那套 FTOOLS 流程）
# ================================

def find_or_make_oot_events(obs_dir: Path) -> Path:
    """
    寻找 OOT 事件文件：
      - 若存在 *OOEVLI*.FIT，则直接使用
      - 否则自动执行 epchain withoutoftime=yes，再搜索一次
    """
    candidates = list(obs_dir.glob("*OOEVLI*.FIT"))
    if candidates:
        if len(candidates) > 1:
            print("[OOT] 警告：发现多个 OOEVLI 文件，将使用：", candidates[0].name)
        else:
            print("[OOT] 自动找到 OOT 事件文件：", candidates[0].name)
        return candidates[0]

    print("[OOT] 未找到 OOEVLI 文件，将在该目录下执行：epchain withoutoftime=yes")
    run_cmd(["epchain", "withoutoftime=yes"], obs_dir)

    candidates = list(obs_dir.glob("*OOEVLI*.FIT"))
    if not candidates:
        raise RuntimeError("[OOT] 运行 epchain 后依然未找到 *OOEVLI*.FIT，请检查 ODF/SAS_ODF 设置。")
    if len(candidates) > 1:
        print("[OOT] 警告：发现多个 OOEVLI 文件，将使用：", candidates[0].name)
    else:
        print("[OOT] 自动找到 OOT 事件文件：", candidates[0].name)
    return candidates[0]


def apply_oot_inplace_with_ftools(
    obs_dir: Path,
    pn_clean: Path,
    src_spec: Path,
    bg_spec: Path,
    xs: float, ys: float, rs: float,
    xb: float, yb: float, rb: float,
    oot_frac: float = 0.063,
):
    """
    按你笔记的做法，对 PNsource_spectrum / PNbackground_spectrum 在原文件内做 OOT 校正：

      1) 找到或生成 OOT 事件文件（OOEVLI）
      2) 从 OOT 事件中提取源/背景 OOT 光谱（区域与 src/bg 一致）
      3) 对 OOT 源/背景光谱做 backscale（可选）
      4) 对每个（源/背景）：
         - 用 fparkey 把 OOT 谱的 COUNTS（TTYPE2）重命名成 CTS_OOT
         - 用 faddcol 把 CTS_OOT 列加到主光谱中
         - 用 fcalc 将 CTS_OOT *= oot_frac
         - 用 fcalc 将 COUNTS = COUNTS - CTS_OOT
      5) 最终：
         - PNsource_spectrum.fits 和 PNbackground_spectrum.fits 本身即为减完 OOT 的光谱
         - 文件名不变，头关键字完整，rmfgen/arfgen 正常使用
    """
    cwd = obs_dir
    oot_events = find_or_make_oot_events(obs_dir)

    oot_src_spec = obs_dir / "PNsource_OOT_spectrum.fits"
    oot_bg_spec  = obs_dir / "PNbackground_OOT_spectrum.fits"

    print("\n[OOT-1] 从 OOT 事件中提取源/背景 OOT 光谱（同一 src/bg 区域）")

    # 源 OOT 光谱
    cmd = [
        "evselect",
        f"table={oot_events.name}",
        "withspectrumset=yes",
        f"spectrumset={oot_src_spec.name}",
        "energycolumn=PI",
        "spectralbinsize=5",
        "withspecranges=yes",
        "specchannelmin=0",
        "specchannelmax=20479",
        f"expression=(FLAG==0)&&(PATTERN<=4)&&((X,Y) IN circle({xs},{ys},{rs}))",
    ]
    run_cmd(cmd, cwd)

    # 背景 OOT 光谱
    cmd = [
        "evselect",
        f"table={oot_events.name}",
        "withspectrumset=yes",
        f"spectrumset={oot_bg_spec.name}",
        "energycolumn=PI",
        "spectralbinsize=5",
        "withspecranges=yes",
        "specchannelmin=0",
        "specchannelmax=20479",
        f"expression=(FLAG==0)&&(PATTERN<=4)&&((X,Y) IN circle({xb},{yb},{rb}))",
    ]
    run_cmd(cmd, cwd)

    print("[OOT-2] backscale OOT 源/背景光谱")
    # 这里可以不指定 badpixlocation，OOT 只是用来近似扣除条纹
    run_cmd(["backscale", f"spectrumset={oot_src_spec.name}"], cwd)
    run_cmd(["backscale", f"spectrumset={oot_bg_spec.name}"], cwd)

    print(f"[OOT-3] 用 FTOOLS 在原始光谱内就地扣除 OOT，oot_frac = {oot_frac:.5f}")

    # === 源光谱：PNsource_spectrum.fits ===
    # (1) OOT 光谱的第 2 列改名为 CTS_OOT
    run_cmd([
        "fparkey",
        "value=CTS_OOT",
        f"fitsfile={oot_src_spec.name}+1",
        "keyword=TTYPE2",
        "add=no"
    ], cwd)

    # (2) 将 CTS_OOT 列加到主源光谱里
    run_cmd([
        "faddcol",
        f"infile={src_spec.name}+1",
        f"colfile={oot_src_spec.name}+1",
        "colname=CTS_OOT"
    ], cwd)

    # (3) 在主源光谱内：CTS_OOT *= oot_frac
    run_cmd([
        "fcalc",
        f"infile={src_spec.name}+1",
        f"outfile={src_spec.name}",
        "clobber=yes",
        "clname=CTS_OOT",
        f"expr=CTS_OOT*{oot_frac}"
    ], cwd)

    # (4) 在主源光谱内：COUNTS = COUNTS - CTS_OOT
    run_cmd([
        "fcalc",
        f"infile={src_spec.name}+1",
        f"outfile={src_spec.name}",
        "clobber=yes",
        "clname=COUNTS",
        "expr=COUNTS-CTS_OOT"
    ], cwd)

    # === 背景光谱：PNbackground_spectrum.fits ===
    run_cmd([
        "fparkey",
        "value=CTS_OOT",
        f"fitsfile={oot_bg_spec.name}+1",
        "keyword=TTYPE2",
        "add=no"
    ], cwd)

    run_cmd([
        "faddcol",
        f"infile={bg_spec.name}+1",
        f"colfile={oot_bg_spec.name}+1",
        "colname=CTS_OOT"
    ], cwd)

    run_cmd([
        "fcalc",
        f"infile={bg_spec.name}+1",
        f"outfile={bg_spec.name}",
        "clobber=yes",
        "clname=CTS_OOT",
        f"expr=CTS_OOT*{oot_frac}"
    ], cwd)

    run_cmd([
        "fcalc",
        f"infile={bg_spec.name}+1",
        f"outfile={bg_spec.name}",
        "clobber=yes",
        "clname=COUNTS",
        "expr=COUNTS-CTS_OOT"
    ], cwd)

    print("[OOT-4] OOT 校正完成：")
    print("         - PNsource_spectrum.fits 已减去 OOT 贡献")
    print("         - PNbackground_spectrum.fits 已减去 OOT 贡献")
    print("         - 文件名不变，后续 RMF/ARF/specgroup 直接使用即可")


# ================================
# 主流程：光谱提取 + OOT + 响应 + 分组
# ================================

def run_module4_for_obsid(obs_dir: Path):
    """
    模块4主体：
      - 设置 SAS 环境变量
      - 检查 PNclean.fits
      - 生成/检查 PNimage.fits
      - 自动调用 ds9，等待 src.reg / bg.reg
      - 解析区域
      - 提取源/背景光谱 + backscale
      - 交互询问是否做 OOT 校正（如是，按你笔记的 FTOOLS 流程就地修改光谱）
      - rmfgen / arfgen
      - specgroup
    """
    print("=" * 60)
    print(f"[模块4] PN 光谱处理：{obs_dir}")
    print("=" * 60)

    # 0) 设置 SAS 环境变量
    setup_sas_env(obs_dir)

    # 1) 检查 PNclean.fits
    pn_clean = obs_dir / "PNclean.fits"
    if not pn_clean.exists():
        raise RuntimeError(f"未找到 PNclean.fits，请先完成模块2+3：{pn_clean}")

    # 2) 图像 + 手动画圈（自动起 ds9）
    src_reg, bg_reg = ensure_image_and_regions(obs_dir, pn_clean)

    print("[4-3] 解析 src.reg")
    xs, ys, rs = parse_ds9_region(src_reg)
    print(f"      src = circle({xs}, {ys}, {rs})")

    print("[4-4] 解析 bg.reg")
    xb, yb, rb = parse_ds9_region(bg_reg)
    print(f"      bg  = circle({xb}, {yb}, {rb})")

    cwd = obs_dir

    # 3) 源光谱（PNclean）
    src_spec = obs_dir / "PNsource_spectrum.fits"
    print("[4-5] 提取源光谱 → PNsource_spectrum.fits")
    cmd = [
        "evselect",
        f"table={pn_clean.name}",
        "withspectrumset=yes",
        f"spectrumset={src_spec.name}",
        "energycolumn=PI",
        "spectralbinsize=5",
        "withspecranges=yes",
        "specchannelmin=0",
        "specchannelmax=20479",
        f"expression=(FLAG==0)&&(PATTERN<=4)&&((X,Y) IN circle({xs},{ys},{rs}))",
    ]
    run_cmd(cmd, cwd)

    # 4) 背景光谱（PNclean）
    bg_spec = obs_dir / "PNbackground_spectrum.fits"
    print("[4-6] 提取背景光谱 → PNbackground_spectrum.fits")
    cmd = [
        "evselect",
        f"table={pn_clean.name}",
        "withspectrumset=yes",
        f"spectrumset={bg_spec.name}",
        "energycolumn=PI",
        "spectralbinsize=5",
        "withspecranges=yes",
        "specchannelmin=0",
        "specchannelmax=20479",
        f"expression=(FLAG==0)&&(PATTERN<=4)&&((X,Y) IN circle({xb},{yb},{rb}))",
    ]
    run_cmd(cmd, cwd)

    # 5) backscale（原始源/背景）
    print("[4-7] backscale (源 + 背景)")
    run_cmd(
        ["backscale", f"spectrumset={src_spec.name}", f"badpixlocation={pn_clean.name}"],
        cwd,
    )
    run_cmd(
        ["backscale", f"spectrumset={bg_spec.name}", f"badpixlocation={pn_clean.name}"],
        cwd,
    )

    # 6) 问你要不要做 OOT
    ans = input("\n❓ 是否对本次 PN 光谱进行 OOT 校正（按你笔记方式）？[y/N]: ").strip().lower()
    apply_oot = (ans == "y")

    if apply_oot:
        print("\n========== [OOT] 已选择进行 OOT 校正（FTOOLS 就地修改 COUNTS） ==========")
        # 这里默认用 Full Frame 的 0.063，你以后可以在这里改数值
        apply_oot_inplace_with_ftools(
            obs_dir=obs_dir,
            pn_clean=pn_clean,
            src_spec=src_spec,
            bg_spec=bg_spec,
            xs=xs, ys=ys, rs=rs,
            xb=xb, yb=yb, rb=rb,
            oot_frac=0.063,
        )
    else:
        print("\n[OOT] 本次不做 OOT 校正，将直接使用未校正的源/背景光谱。")

    # 7) RMF（用当前的 PNsource_spectrum.fits，无论是否已减 OOT）
    rmf_file = obs_dir / "PN.rmf"
    print("[4-8] 生成 RMF → PN.rmf")
    cmd = ["rmfgen", f"spectrumset={src_spec.name}", f"rmfset={rmf_file.name}"]
    run_cmd(cmd, cwd)

    # 8) ARF
    arf_file = obs_dir / "PN.arf"
    print("[4-9] 生成 ARF → PN.arf")
    cmd = [
        "arfgen",
        f"spectrumset={src_spec.name}",
        f"arfset={arf_file.name}",
        "withrmfset=yes",
        f"rmfset={rmf_file.name}",
        f"badpixlocation={pn_clean.name}",
        "detmaptype=psf",
        "applyabsfluxcorr=yes",
    ]
    run_cmd(cmd, cwd)

    # 9) 分组（直接用当前 PNsource_spectrum / PNbackground_spectrum）
    grouped = obs_dir / "PN_spectrum_grp.fits"
    print("[4-10] specgroup → PN_spectrum_grp.fits")
    cmd = [
        "specgroup",
        f"spectrumset={src_spec.name}",
        "mincounts=25",
        "oversample=3",
        f"rmfset={rmf_file.name}",
        f"arfset={arf_file.name}",
        f"backgndset={bg_spec.name}",
        f"groupedset={grouped.name}",
    ]
    run_cmd(cmd, cwd)

    print("\n==========================")
    print("模块4 完成！光谱文件已生成：")
    print(f"  源光谱：   {src_spec}")
    print(f"  背景光谱： {bg_spec}")
    print(f"  RMF：      {rmf_file}")
    print(f"  ARF：      {arf_file}")
    print(f"  分组谱：   {grouped}")
    if apply_oot:
        print("  （本次已按笔记方式进行 OOT 校正：COUNTS 已扣除 OOT）")
    else:
        print("  （本次未进行 OOT 校正）")
    print("==========================\n")


# ================================
# CLI — 一次可处理多个 OBSID
# ================================

def main():
    parser = argparse.ArgumentParser(
        description="模块4：PNclean.fits → 图像 + 手动画圈 + 光谱提取 + （交互式可选 OOT 校正）+ RMF/ARF + 分组"
    )
    parser.add_argument("source", help="源名，例如 IRASF12397+3333")
    parser.add_argument(
        "obsids",
        nargs="+",
        help="观测号（一个或多个），例如：0202180201",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "data" / "XMM",
        help="数据根路径（默认 ~/data/XMM）",
    )

    args = parser.parse_args()
    data_root = args.data_root.expanduser()

    for obsid in args.obsids:
        obs_dir = data_root / args.source / obsid
        if not obs_dir.exists():
            print(f"[错误] 目录不存在，跳过：{obs_dir}", file=sys.stderr)
            continue

        try:
            run_module4_for_obsid(obs_dir)
        except Exception as e:
            print(f"[错误] 处理 {obsid} 时失败：{e}", file=sys.stderr)


if __name__ == "__main__":
    main()
