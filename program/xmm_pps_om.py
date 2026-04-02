#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import shutil
import tarfile
import requests

# =========================
# 配置
# =========================
BASE_URL = "https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio"
OM_RSP_DIR = os.path.expanduser("~/data/XMM/om_effarea_v2.0")

CANDIDATE_PRODUCTS = [
    "OMCOMBOBSMLI0000",
    "OMCOMBOBSMER0000",
    "OMX000OBSMLI0000",
]


# =========================
# 基础工具
# =========================
def build_url(obsid: str, stem: str) -> str:
    filename = f"P{obsid}{stem}.FTZ"
    return (
        f"{BASE_URL}?obsno={obsid}"
        f"&extension=FTZ"
        f"&name={filename}"
    )


def maybe_html(filepath: str) -> bool:
    with open(filepath, "rb") as f:
        head = f.read(512).lower()
    return (b"<html" in head) or (b"<!doctype html" in head)


def looks_like_gzip(filepath: str) -> bool:
    with open(filepath, "rb") as f:
        sig = f.read(2)
    return sig == b"\x1f\x8b"


def is_fits_file(filepath: str) -> bool:
    try:
        with open(filepath, "rb") as f:
            head = f.read(80)
        return head.startswith(b"SIMPLE  ") or head.startswith(b"XTENSION")
    except Exception:
        return False


def download_one(url: str, output_file: str, timeout: int = 120) -> bool:
    try:
        r = requests.get(url, stream=True, timeout=timeout, allow_redirects=True)
        r.raise_for_status()
    except requests.RequestException as e:
        print(f"[ERROR] Request failed: {e}")
        return False

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)

    size = os.path.getsize(output_file)
    print(f"[INFO] Saved: {output_file} ({size / 1024:.1f} KB)")

    if size == 0:
        print("[WARN] Empty file.")
        return False

    if maybe_html(output_file):
        print("[WARN] Downloaded content looks like HTML, not FTZ/FITS.")
        return False

    if not looks_like_gzip(output_file):
        print("[WARN] File is not gzip-compressed.")
        return False

    return True


def gunzip_to_file(ftz_path: str, out_path: str) -> None:
    with gzip.open(ftz_path, "rb") as fin, open(out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    print(f"[INFO] Decompressed to: {out_path}")


def extract_tar(tar_path: str, extract_dir: str) -> None:
    with tarfile.open(tar_path, "r") as tar:
        tar.extractall(path=extract_dir)
    print(f"[INFO] Extracted tar archive: {tar_path} -> {extract_dir}")


def copy_response_files(target_dir: str) -> None:
    if not os.path.isdir(OM_RSP_DIR):
        print(f"[WARN] Response dir not found: {OM_RSP_DIR}")
        return

    copied = 0
    for fname in os.listdir(OM_RSP_DIR):
        if fname.lower().endswith((".rsp", ".rmf")):
            src = os.path.join(OM_RSP_DIR, fname)
            dst = os.path.join(target_dir, fname)
            shutil.copy2(src, dst)
            copied += 1

    print(f"[INFO] Copied {copied} response files.")


# =========================
# 搜索真正的 OM source list
# =========================
def score_candidate(name: str, obsid: str) -> int:
    """
    分数越高越优先
    """
    up = name.upper()
    score = 0

    if f"P{obsid}OMCOMBOBSMLI0000.FIT".upper() == up:
        score += 100
    elif f"P{obsid}OMCOMBOBSMER0000.FIT".upper() == up:
        score += 90
    elif f"P{obsid}OMX000OBSMLI0000.FIT".upper() == up:
        score += 80

    if "OMCOMBOBSMLI" in up:
        score += 50
    elif "OMCOMBOBSMER" in up:
        score += 40
    elif "OBSMLI" in up:
        score += 30

    if up.endswith(".FIT") or up.endswith(".FITS"):
        score += 10

    return score


def find_real_srclist(search_dir: str, obsid: str) -> str | None:
    candidates = []

    for root, _, files in os.walk(search_dir):
        for name in files:
            up = name.upper()
            if not up.endswith((".FIT", ".FITS")):
                continue
            if not any(k in up for k in ["OMCOMBOBSMLI", "OMCOMBOBSMER", "OBSMLI"]):
                continue

            full = os.path.join(root, name)
            if is_fits_file(full):
                candidates.append((score_candidate(name, obsid), full))

    if not candidates:
        return None

    candidates.sort(reverse=True, key=lambda x: x[0])
    return candidates[0][1]


def extract_nested_tars(search_dir: str, max_rounds: int = 3) -> None:
    """
    自动多轮解 tar，防止 archive 里还有 tar
    """
    for _ in range(max_rounds):
        found_tar = False
        for root, _, files in os.walk(search_dir):
            for name in files:
                full = os.path.join(root, name)
                try:
                    if tarfile.is_tarfile(full):
                        outdir = full + "_dir"
                        if not os.path.exists(outdir):
                            os.makedirs(outdir, exist_ok=True)
                            extract_tar(full, outdir)
                        found_tar = True
                except Exception:
                    pass
        if not found_tar:
            break


# =========================
# 清理中间文件
# =========================
def cleanup_intermediate_files(target_dir: str, keep_files: set[str]) -> None:
    for entry in os.listdir(target_dir):
        full = os.path.join(target_dir, entry)

        if full in keep_files:
            continue

        # 保留 response
        if entry.lower().endswith((".rsp", ".rmf", ".pha")):
            continue

        # 删除中间目录
        if os.path.isdir(full) and (
            entry.endswith("_extracted")
            or entry.endswith("_dir")
        ):
            shutil.rmtree(full, ignore_errors=True)
            print(f"[CLEAN] Removed dir: {full}")
            continue

        # 删除中间文件
        if os.path.isfile(full) and (
            entry.endswith(".FTZ")
            or entry.endswith(".unpacked")
            or entry.endswith(".TAR")
        ):
            os.remove(full)
            print(f"[CLEAN] Removed file: {full}")


# =========================
# 主流程
# =========================
def download_and_prepare_om(obsid: str, target_dir: str, copy_rsp: bool = True) -> str | None:
    os.makedirs(target_dir, exist_ok=True)

    final_srclist = None
    keep_files = set()

    for stem in CANDIDATE_PRODUCTS:
        ftz_name = f"P{obsid}{stem}.FTZ"
        ftz_path = os.path.join(target_dir, ftz_name)
        unpacked_path = os.path.join(target_dir, f"{stem}.unpacked")
        extract_dir = os.path.join(target_dir, f"{stem}_extracted")
        url = build_url(obsid, stem)

        print("\n" + "=" * 70)
        print(f"[INFO] Trying product: {stem}")
        print(f"[INFO] URL: {url}")

        ok = download_one(url, ftz_path)
        if not ok:
            if os.path.exists(ftz_path):
                os.remove(ftz_path)
            continue

        gunzip_to_file(ftz_path, unpacked_path)

        # 直接是 FITS
        if is_fits_file(unpacked_path):
            final_srclist = os.path.join(target_dir, f"P{obsid}{stem}.FIT")
            os.replace(unpacked_path, final_srclist)
            keep_files.add(final_srclist)
            print(f"[OK] Direct FITS source list found: {final_srclist}")
            break

        # 是 tar 包
        if tarfile.is_tarfile(unpacked_path):
            os.makedirs(extract_dir, exist_ok=True)
            extract_tar(unpacked_path, extract_dir)
            extract_nested_tars(extract_dir, max_rounds=3)

            found = find_real_srclist(extract_dir, obsid)
            if found is not None:
                final_srclist = os.path.join(target_dir, os.path.basename(found))
                if os.path.abspath(found) != os.path.abspath(final_srclist):
                    shutil.copy2(found, final_srclist)
                keep_files.add(final_srclist)
                print(f"[OK] Real source list found: {final_srclist}")
                break
            else:
                print(f"[WARN] No usable source list found inside {stem}")

    if final_srclist is None:
        print("[ERROR] No suitable OM source list product was found.")
        return None

    if copy_rsp:
        copy_response_files(target_dir)

    cleanup_intermediate_files(target_dir, keep_files)

    return final_srclist


def print_om2pha_command(srclist_fit: str, ra: float, dec: float, output_pha: str) -> None:
    print("\nUse this om2pha command:")
    print(
        f"om2pha srclist={os.path.basename(srclist_fit)} "
        f"ra={ra:.6f} dec={dec:.6f} output={output_pha}"
    )


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 xmm_pps_om.py <source_name> <obsid>")
        sys.exit(1)

    source_name = sys.argv[1]
    obsid = sys.argv[2]

    coords = {
        "RE1034+396": (180.27824, 59.06459),
    }

    if source_name not in coords:
        print(f"[ERROR] No coordinates configured for {source_name}")
        sys.exit(1)

    ra, dec = coords[source_name]
    target_dir = os.path.expanduser(f"~/data/XMM/{source_name}/{obsid}")
    output_pha = f"{source_name.replace('+', '_')}_{obsid}_om.pha"

    srclist_fit = download_and_prepare_om(obsid, target_dir, copy_rsp=True)
    if srclist_fit is not None:
        print_om2pha_command(srclist_fit, ra, dec, output_pha)