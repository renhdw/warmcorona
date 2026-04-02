#!/usr/bin/env python3
"""
模块 1（代理模式）：使用代理下载 XMM ODF （仅 ODF，不含 PPS）

功能：
  对给定源名和多个 OBSID：
    - 用 NXSA &level=ODF 只下载 ODF
    - 使用当前 shell 配置的代理（http_proxy / https_proxy）
    - wget 断点续传 + 超时 + 重试
    - 若本地已存在非零 obsid.tar，则跳过下载，直接解压
    - 解压 obsid.tar 前后比较目录，识别新解出的内容
    - 统一整理到 <SOURCE>/<OBSID>/ 目录
    - 若 ODF 内还有 *.TAR，再在 <OBSID>/ 目录内继续解压并删除内层 TAR
    - 删除外层 obsid.tar 节省空间
"""

import argparse
import os
import subprocess
from pathlib import Path


def run_cmd(cmd, cwd: Path, env=None):
    """在 cwd 下运行 cmd，失败则抛异常。"""
    print(f"[CMD@{cwd.name}] {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=str(cwd), env=env)
    if r.returncode != 0:
        raise RuntimeError(f"命令失败: {' '.join(cmd)} (returncode={r.returncode})")


def download_odf_proxy(source_dir: Path, obsid: str, proxy_env: dict):
    """
    使用代理下载单个 obsid 的 ODF 包，并进行两层解压。

    最终结构期望：
        <source_dir>/<obsid>/
            SUM.SAS
            *0000.FIT
            其他 ODF 文件
    """
    print("=" * 60)
    print(f"源目录: {source_dir} | OBSID: {obsid}")
    print("=" * 60)

    source_dir.mkdir(parents=True, exist_ok=True)

    # 只下载 ODF
    url = (
        "https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio"
        f"?obsno={obsid}&level=ODF"
    )
    out_tar = source_dir / f"{obsid}.tar"

    # ==== 若本地已有非零 tar，直接解压，不再走代理 ====
    if out_tar.exists() and out_tar.stat().st_size > 0:
        print(f"[{obsid}] 检测到本地已存在非零 tar：{out_tar}")
        print(f"        将跳过下载，直接使用该文件进行解压。")
    else:
        # wget 通过代理，下大文件尽量稳一点
        wget_cmd = [
            "wget",
            "-c",              # 断点续传
            "--tries=5",       # 重试 5 次
            "--timeout=300",   # 单次连接超时 300 s
            "-O", str(out_tar),
            url,
        ]

        print(f"[{obsid}] 开始下载：{url}")
        run_cmd(wget_cmd, cwd=source_dir, env=proxy_env)

        if not out_tar.exists() or out_tar.stat().st_size == 0:
            raise RuntimeError(f"[{obsid}] 下载的 {out_tar} 不存在或大小为 0")

        print(f"[{obsid}] 下载完成：{out_tar}")

    # ==== 解压前记录一下已有的内容，用于识别新解压出来的东西 ====
    before_items = {p.name for p in source_dir.iterdir()}

    # 第一层解压：在 source_dir 下解开 obsid.tar
    print(f"解压 {out_tar.name} ...")
    run_cmd(["tar", "xf", out_tar.name], cwd=source_dir, env=proxy_env)

    # 解压后目录内容
    after_items = {p.name for p in source_dir.iterdir()}
    new_items = after_items - before_items  # 这就是 tar 新解出的东西

    # 目标 obsid 目录
    obs_dir = source_dir / obsid
    if not obs_dir.exists():
        obs_dir.mkdir()
        print(f"创建观测目录：{obs_dir}")

    # 1) 先把新解出来的目录搬进 obsid 目录
    for name in sorted(new_items):
        path = source_dir / name
        if path.is_dir() and name != obsid:
            print(f"  移动目录 {name} → {obs_dir.name}/")
            path.rename(obs_dir / name)

    # 2) 再把新解出来的文件搬进 obsid 目录
    for name in sorted(new_items):
        path = source_dir / name
        if not path.exists():
            continue
        if path.is_file():
            print(f"  移动文件 {name} → {obs_dir.name}/")
            path.rename(obs_dir / name)

    print(f"✓ 观测目录内容已归档至：{obs_dir}")

    # 第二层：在 obs_dir 下解压可能存在的 *.TAR（ODF 内层包）
    inner_tars = list(obs_dir.glob("*.TAR"))
    for inner in inner_tars:
        print(f"解压内层 {inner.name} ...")
        run_cmd(["tar", "xf", inner.name], cwd=obs_dir, env=proxy_env)
        print(f"删除内层 {inner.name}")
        inner.unlink()

    # 删除外层 tar
    print(f"删除外层 {out_tar.name}")
    out_tar.unlink(missing_ok=True)

    print(f"[{obsid}] 处理完毕。\n")


def main():
    parser = argparse.ArgumentParser(
        description="使用代理从 NXSA 下载 XMM ODF（仅 ODF），并解压到 ~/data/XMM/<SOURCE>/<OBSID>/"
    )
    parser.add_argument("source", help="源名，例如 1H0419-577")
    parser.add_argument("obsids", nargs="+", help="观测号列表")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "data" / "XMM",
        help="数据根目录（默认 ~/data/XMM）",
    )
    args = parser.parse_args()

    data_root = args.data_root.expanduser()
    source = args.source
    obsids = args.obsids

    print(f"数据根目录: {data_root}")
    print(f"源名字: {source}")
    print(f"观测号列表: {', '.join(obsids)}")

    # 继承当前 shell 的代理设置
    proxy_env = os.environ.copy()
    if "http_proxy" not in proxy_env and "https_proxy" not in proxy_env:
        print("⚠️ 当前环境没有设置 http_proxy/https_proxy，仍将继续尝试，但可能不会走代理。")

    source_dir = data_root / source

    for obsid in obsids:
        try:
            download_odf_proxy(source_dir, obsid, proxy_env)
        except Exception as e:
            print(f"❌ 观测 {obsid} 处理失败：{e}\n    将继续处理后续观测。\n")


if __name__ == "__main__":
    main()
