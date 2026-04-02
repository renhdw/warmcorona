#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import urllib.parse as up
import re

from astroquery.simbad import Simbad
import requests


# =========================
# 1. SIMBAD 查询部分
# =========================

def query_simbad_basic(name: str):
    """
    查询 SIMBAD：
      - main_id
      - ra, dec (deg)
      - 红移 rvz_redshift
    """

    custom = Simbad()
    custom.add_votable_fields("rvz_redshift")

    result = custom.query_object(name)
    if result is None or len(result) == 0:
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
    构造 HEASARC w3nh 的参数字典

    按你提供的 URL 一模一样：
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
        "raw_text": str                  # 原始返回内容（调试用）
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
# 3. 主程序
# =========================

def main():
    if len(sys.argv) < 2:
        print("用法: python3 test_simbad_nh.py '1H0419-577'")
        sys.exit(1)

    src_name = sys.argv[1]

    # ---- SIMBAD ----
    print(f"🔍 查询 SIMBAD: {src_name} ...")
    info = query_simbad_basic(src_name)

    print("\n=== SIMBAD 结果 ===")
    print("  MAIN_ID:", info["main_id"])
    print(f"  RA  (deg): {info['ra_deg']:.6f}")
    print(f"  Dec (deg): {info['dec_deg']:.6f}")
    if info["z"] is None:
        print("  Redshift: 无（SIMBAD 未给出可用红移）")
    else:
        print(f"  Redshift: {info['z']:.6f}")

    print("\n  SIMBAD 返回列名：")
    print("   ", info["raw_colnames"])

    print("\n🔗 SIMBAD 页面：")
    print("   ", build_simbad_url(src_name))

    # ---- w3nh ----
    print("\n🔍 查询 w3nh 平均 N_H ...")
    nh_info = fetch_w3nh_nh(info["ra_deg"], info["dec_deg"], radius_deg=0.1)

    avg_nh = nh_info["avg_nh"]
    wavg_nh = nh_info["wavg_nh"]

    print("\n=== w3nh 结果 ===")
    if avg_nh is None and wavg_nh is None:
        print("  未能从 w3nh 页面解析出 Average nH / Weighted average nH")
        # 调试：你如果想看原始内容，可以手动写入一个文件
        # with open("w3nh_debug.html", "w", encoding="utf-8") as f:
        #     f.write(nh_info["raw_text"])
    else:
        if avg_nh is not None:
            print(f"  Average nH (cm^-2): {avg_nh:.3e}")
            print(f"  → TBabs nH = {avg_nh/1e22:.3e} (10^22 cm^-2)")
        if wavg_nh is not None:
            print(f"  Weighted average nH (cm^-2): {wavg_nh:.3e}")
            print(f"  → TBabs nH = {wavg_nh/1e22:.3e} (10^22 cm^-2)")

    print("\n🔗 w3nh 页面（可人工检查）：")
    print("   ", nh_info["url"])


if __name__ == "__main__":
    main()
