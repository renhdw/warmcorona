#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模块 5：自动查询 z / N_H 并生成 XSPEC xcm 文件（localmodel 版，纯净版）

功能：
  - 输入：源名（用于 SIMBAD / 路径） + 若干 OBSID
  - 自动：
      1) 用 SIMBAD 查询源的 redshift z、RA/Dec（带名字变体尝试）
      2) 用 HEASARC w3nh 根据 RA/Dec 查询 Galactic N_H
      3) 为每个 OBSID 生成 slab / sphere 两个 xcm 文件
         * xcm 内容与原来模板一致（不含额外网页注释）
         * 使用本地模型 warmcomslab / warmcomsphere
         * 自动设置参数列表中的：
             p1 = TBabs:nH = N_H / 1e22
             p3 = zTBabs:Redshift = z

目录结构约定：
  ~/data/XMM/<SOURCE>/<OBSID>/
      PN_spectrum_grp.fits   # 模块4生成
      xcm/
          fit_model_Tom_slab_<SOURCE>_<OBSID>.xcm
          fit_model_tv_Tom_sphere_<SOURCE>_<OBSID>.xcm

另外在：
  ~/data/XMM/<SOURCE>/xcm_info_<SOURCE>.log
  记录 SIMBAD / w3nh 查询结果（方便以后查 z / N_H / URL）

使用示例：
  python3 xmm_make_xcm_localmodel.py 1H0419-577 0604720301 0148000401
"""

import sys
import re
import urllib.parse as up
from pathlib import Path

import requests
from astroquery.simbad import Simbad


# =========================
# 1. SIMBAD 查询部分（带名字变体）
# =========================

def generate_name_candidates(name: str):
    """
    根据 XMM 目标名，生成一组可能的 SIMBAD 名称候选。
    例如：
      IRASF12397+3333  -> ["IRASF12397+3333", "IRAS F12397+3333", "IRAS 12397+3333"]
    后面如果遇到其他 pattern，可以继续往这里加规则。
    """
    cand = []
    n0 = name.strip()
    cand.append(n0)

    upper = n0.upper().replace(" ", "")

    # 针对 IRASFxxxx 的特殊处理
    if upper.startswith("IRASF"):
        # 原样 IRASF12397+3333
        rest = n0.strip()[5:]  # 去掉 "IRASF"
        # 变体 1: IRAS F12397+3333
        cand.append(f"IRAS F{rest}")
        # 变体 2: IRAS 12397+3333
        cand.append(f"IRAS {rest}")

    # 针对 IRASxxxx（中间没空格）的处理
    elif upper.startswith("IRAS") and " " not in n0[:5]:
        rest = n0.strip()[4:]  # 去掉 "IRAS"
        cand.append(f"IRAS {rest}")

    # 可以在这里继续添加其他规则，例如：
    # if upper.startswith("1H"):
    #     ...

    # 去重保持顺序
    seen = set()
    out = []
    for c in cand:
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out


def query_simbad_basic(name: str):
    """
    查询 SIMBAD：
      - main_id
      - ra, dec (deg)
      - 红移 rvz_redshift

    会自动尝试常见名称变体（比如 IRASFxxxx → IRAS Fxxxx）。
    返回字典包含：
      main_id, ra_deg, dec_deg, z, raw_colnames, used_name
    """
    custom = Simbad()
    custom.add_votable_fields("rvz_redshift")

    last_err = None
    result = None
    used_name = None

    for cand in generate_name_candidates(name):
        try:
            result = custom.query_object(cand)
        except Exception as e:
            last_err = e
            result = None

        if result is not None and len(result) > 0:
            used_name = cand
            break

    if result is None or len(result) == 0:
        if last_err is not None:
            raise RuntimeError(f"SIMBAD 中没有找到对象：{name}，最后尝试错误：{last_err}")
        raise RuntimeError(f"SIMBAD 中没有找到对象：{name}")

    row = result[0]

    # main_id
    if "main_id" in row.colnames:
        mid = row["main_id"]
        if hasattr(mid, "decode"):
            mid = mid.decode("utf-8")
        main_id = str(mid)
    else:
        main_id = name

    # ra / dec：已经是度
    try:
        ra_deg = float(row["ra"])
        dec_deg = float(row["dec"])
    except Exception as e:
        raise RuntimeError(f"无法从 SIMBAD 结果读取 ra/dec: {e}")

    # 红移：rvz_redshift
    z_val = None
    if "rvz_redshift" in row.colnames:
        zv = row["rvz_redshift"]
        if zv is not None:
            try:
                z_val = float(zv)
            except Exception:
                z_val = None

    return {
        "main_id": main_id,
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "z": z_val,
        "raw_colnames": row.colnames,
        "used_name": used_name or name,
    }


def build_simbad_url(name: str) -> str:
    base = "https://simbad.cds.unistra.fr/simbad/sim-basic"
    q = up.urlencode({"Ident": name})
    return f"{base}?{q}"


# =========================
# 2. w3nh 查询 + 解析部分
# =========================

def build_w3nh_params(ra_deg: float, dec_deg: float, radius_deg: float = 0.1) -> dict:
    """
    按验证过的 URL 构造 w3nh 查询参数：
    https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl
      ?Entry=66.50300,+-57.20049
      &NR=GRB/SIMBAD+Sesame/NED
      &CoordSys=Equatorial
      &equinox=2000
      &radius=0.1
      &usemap=0
    """
    entry = f"{ra_deg:.5f}, {dec_deg:.5f}"

    params = {
        "Entry": entry,
        "NR": "GRB/SIMBAD+Sesame/NED",
        "CoordSys": "Equatorial",
        "equinox": "2000",
        "radius": f"{radius_deg}",
        "usemap": "0",
    }
    return params


def build_w3nh_url_from_params(params: dict) -> str:
    base = "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl"
    return f"{base}?{up.urlencode(params)}"


_float_re = re.compile(r"([0-9.+-]+E[+-]?[0-9]+)", re.IGNORECASE)


def _extract_scientific_number(line: str):
    """
    在一行里通过正则提取第一个科学计数法数字，如 1.14E+20。
    没有则返回 None。
    """
    m = _float_re.search(line)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def fetch_w3nh_nh(ra_deg: float, dec_deg: float, radius_deg: float = 0.1):
    """
    调用 w3nh CGI，解析文本中的 Average nH / Weighted average nH。

    返回：
      {
        "avg_nh": float or None,         # cm^-2
        "wavg_nh": float or None,        # cm^-2
        "url": str,                      # 对应这次查询的网页 URL
        "raw_text": str                  # 原始返回内容（可选调试）
      }
    """
    base = "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl"
    params = build_w3nh_params(ra_deg, dec_deg, radius_deg=radius_deg)
    url = build_w3nh_url_from_params(params)

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0 Safari/537.36"
        )
    }

    resp = requests.get(base, params=params, timeout=20, headers=headers)
    resp.raise_for_status()

    text = resp.text

    avg_nh = None
    wavg_nh = None

    for line in text.splitlines():
        low = line.lower()

        # Average nH 行（不含 weighted）
        if "average nh" in low and "weighted" not in low:
            val = _extract_scientific_number(line)
            if val is not None:
                avg_nh = val

        # Weighted average nH 行
        if "weighted average nh" in low:
            val = _extract_scientific_number(line)
            if val is not None:
                wavg_nh = val

    return {
        "avg_nh": avg_nh,
        "wavg_nh": wavg_nh,
        "url": url,
        "raw_text": text,
    }


# =========================
# 3. XSPEC xcm 生成部分（localmodel）
# =========================

def build_params_for_source(z: float | None, nh_gal_1e22: float | None):
    """
    根据 redshift 和 Galactic N_H (10^22 cm^-2) 构造参数列表。
    其余参数用你之前的一组默认初值。

    当前模型顺序为：
      model TBabs*zTBabs(nthComp + warmcomXXX + xillver)

    对应参数顺序：
      1  TBabs:nH
      2  zTBabs:nH
      3  zTBabs:Redshift
      4  nthComp:Gamma
      5  nthComp:kT_e
      6  nthComp:kT_bb
      7  nthComp:inp_type
      8  nthComp:Redshift
      9  nthComp:norm
      10 warmcom:te
      11 warmcom:tau
      12 warmcom:Redshift
      13 warmcom:norm
      14 xillver:gamma
      15 xillver:Afe
      16 xillver:Ecut
      17 xillver:logxi
      18 xillver:z
      19 xillver:Incl
      20 xillver:refl_frac
      21 xillver:norm
    """
    if z is None:
        z = 0.0
    if nh_gal_1e22 is None:
        nh_gal_1e22 = 0.01  # 兜底

    params = [
        nh_gal_1e22,  # 1  TBabs:nH (Galactic N_H)
        1.0,          # 2  zTBabs:nH
        z,            # 3  zTBabs:Redshift
        1.5,          # 4  nthComp:Gamma
        100.0,        # 5  nthComp:kT_e
        3.0e-3,       # 6  nthComp:kT_bb
        0,            # 7  nthComp:inp_type
        z,            # 8  nthComp:Redshift (= p3)
        1.0,          # 9  nthComp:norm
        0.5,          # 10 warmcom:te
        12.0,         # 11 warmcom:tau
        z,            # 12 warmcom:Redshift (= p3)
        1.0,          # 13 warmcom:norm
        1.5,          # 14 xillver:gamma (= nthComp:Gamma)
        1.0,          # 15 xillver:Afe
        300.0,        # 16 xillver:Ecut
        0.0,          # 17 xillver:logxi
        z,            # 18 xillver:z (= p3)
        30.0,         # 19 xillver:Incl
        -1.0,         # 20 xillver:refl_frac
        1.0,          # 21 xillver:norm
    ]
    return params


def generate_xspec_xcm(
    output_filename: str,
    warm_model_name: str,
    parameters: list[float],
    output_dir: Path,
    additional_commands: list[str] | None = None,
):
    """
    生成一个 XSPEC .xcm 文件（使用本地 warmcom 模型）
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_filename

    param_lines = "\n".join(str(p) for p in parameters)
    extra_cmds = "\n".join(additional_commands) if additional_commands else ""

    xcm_content = f"""
data PN_spectrum_grp.fits
cpd /xs
setpl energy
setpl add
ignore 1.8-2.2
ignore **-0.3
ignore 10.0-**
lmod relxill ~/data/monk/plot/warmcorona/model/relxill/relxill_model_v2.3
lmod warmcom ~/data/monk/plot/warmcorona/model/warmcom
model TBabs*zTBabs(nthComp + {warm_model_name} + xillver)
{param_lines}
{extra_cmds}
""".strip() + "\n"

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(xcm_content)

    print(f"✅ XSPEC xcm 文件已生成: {output_path}")


# =========================
# 4. 配置区（localmodel 名称等）
# =========================

WARM_MODEL_SLAB = "warmcomslab"
WARM_MODEL_SPHERE = "warmcomsphere"

XCM_SLAB_TEMPLATE = "fit_model_Tom_slab_{source}_{obsid}.xcm"
XCM_SPHERE_TEMPLATE = "fit_model_tv_Tom_sphere_{source}_{obsid}.xcm"

EXTRA_CMDS = [
    "plot ldata emo re",
    "freeze 1",
    "freeze 5",
    "new 8=3",
    "new 12=3",
    "new 14=4",
    "new 18=3",
    "freeze 15",
    "freeze 16",
    "freeze 17",
    "new 0",
]


# =========================
# 5. 主程序
# =========================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="模块5：自动查询 z / N_H 并生成 XSPEC xcm 文件（slab + sphere，localmodel 版）"
    )
    parser.add_argument("source", help="源名，例如 1H0419-577 （用于 SIMBAD & 路径）")
    parser.add_argument("obsids", nargs="+", help="观测号 OBSID（一个或多个）")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "data" / "XMM",
        help="数据根目录（默认：~/data/XMM）",
    )
    args = parser.parse_args()

    data_root: Path = args.data_root.expanduser()
    source_name: str = args.source
    obsids: list[str] = args.obsids

    print(f"数据根目录: {data_root}")
    print(f"源名字: {source_name}")
    print(f"观测号列表: {', '.join(obsids)}")

    # ---- 1) 查 SIMBAD + w3nh（整源共用一次） ----
    print(f"\n🔍 查询 SIMBAD: {source_name} ...")
    simbad_info = query_simbad_basic(source_name)

    main_id = simbad_info["main_id"]
    ra_deg = simbad_info["ra_deg"]
    dec_deg = simbad_info["dec_deg"]
    z = simbad_info["z"]
    used_name = simbad_info.get("used_name", source_name)

    simbad_url = build_simbad_url(used_name)

    print("\n=== SIMBAD 结果 ===")
    print("  查询名  :", used_name)
    print("  MAIN_ID :", main_id)
    print(f"  RA  (deg): {ra_deg:.6f}")
    print(f"  Dec (deg): {dec_deg:.6f}")
    if z is None:
        print("  Redshift: 无（SIMBAD 未给出可用红移，将使用 z=0.0）")
    else:
        print(f"  Redshift: {z:.6f}")

    print("\n🔍 查询 w3nh 平均 N_H ...")
    nh_info = fetch_w3nh_nh(ra_deg, dec_deg, radius_deg=0.1)
    avg_nh = nh_info["avg_nh"]
    wavg_nh = nh_info["wavg_nh"]
    w3nh_url = nh_info["url"]

    nh_used_cm2 = None

    print("\n=== w3nh 结果 ===")
    if avg_nh is None and wavg_nh is None:
        print("  未能解析 Average nH / Weighted average nH，将使用默认 N_H = 1e20 cm^-2")
        nh_used_cm2 = 1.0e20
    else:
        if avg_nh is not None:
            print(f"  Average nH (cm^-2): {avg_nh:.3e}")
        if wavg_nh is not None:
            print(f"  Weighted average nH (cm^-2): {wavg_nh:.3e}")
        nh_used_cm2 = wavg_nh if wavg_nh is not None else avg_nh

    nh_used_1e22 = nh_used_cm2 / 1e22 if nh_used_cm2 is not None else 0.01
    print(f"\n  用于 TBabs 的 nH = {nh_used_1e22:.3e} (10^22 cm^-2)")
    print("\n🔗 SIMBAD:", simbad_url)
    print("🔗 w3nh  :", w3nh_url)
    print()

    # ---- 在源目录下记一份 log（不写进 xcm） ----
    source_dir = data_root / source_name
    source_dir.mkdir(parents=True, exist_ok=True)
    log_path = source_dir / f"xcm_info_{source_name}.log"
    with open(log_path, "a", encoding="utf-8") as f:
        f.write("\n" + "=" * 60 + "\n")
        f.write(f"Source: {source_name} (SIMBAD MAIN_ID = {main_id})\n")
        f.write(f"Used name : {used_name}\n")
        f.write(f"RA (deg)  = {ra_deg:.6f}\n")
        f.write(f"Dec (deg) = {dec_deg:.6f}\n")
        f.write(f"z         = {0.0 if z is None else z:.6f}\n")
        f.write(f"nH_used   = {nh_used_1e22:.6e} (10^22 cm^-2)\n")
        f.write(f"SIMBAD URL: {simbad_url}\n")
        f.write(f"w3nh URL  : {w3nh_url}\n")
        f.write("ObsIDs: " + ", ".join(obsids) + "\n")

    # ---- 2) 构造参数列表（整源共用） ----
    base_params = build_params_for_source(z, nh_used_1e22)

    # ---- 3) 针对每个 OBSID 生成 xcm ----
    for obsid in obsids:
        obs_dir = data_root / source_name / obsid
        if not obs_dir.exists():
            print(f"[警告] 观测目录不存在，跳过：{obsid} ({obs_dir})", file=sys.stderr)
            continue

        spec_path = obs_dir / "PN_spectrum_grp.fits"
        if not spec_path.exists():
            print(f"[警告] 未找到 PN_spectrum_grp.fits，跳过：{spec_path}", file=sys.stderr)
            continue

        xcm_dir = obs_dir / "xcm"
        print("=" * 60)
        print(f"[模块5] 为 {source_name} / {obsid} 生成 xcm 文件")
        print(f"  观测目录: {obs_dir}")
        print(f"  xcm 目录: {xcm_dir}")
        print("=" * 60)

        params_for_this = list(base_params)

        # slab
        xcm_slab_name = XCM_SLAB_TEMPLATE.format(source=source_name, obsid=obsid)
        generate_xspec_xcm(
            output_filename=xcm_slab_name,
            warm_model_name=WARM_MODEL_SLAB,
            parameters=params_for_this,
            output_dir=xcm_dir,
            additional_commands=EXTRA_CMDS,
        )

        # sphere
        xcm_sphere_name = XCM_SPHERE_TEMPLATE.format(source=source_name, obsid=obsid)
        generate_xspec_xcm(
            output_filename=xcm_sphere_name,
            warm_model_name=WARM_MODEL_SPHERE,
            parameters=params_for_this,
            output_dir=xcm_dir,
            additional_commands=EXTRA_CMDS,
        )

    print("\n[模块5] 全部处理结束。")
    print(f"  查询信息已记录在: {log_path}")


if __name__ == "__main__":
    main()