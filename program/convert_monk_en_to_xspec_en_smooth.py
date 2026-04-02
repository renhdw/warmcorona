#!/usr/bin/env python
# coding: utf-8

# In[1]:


import struct
import os
import numpy as np
import matplotlib.pyplot as plt


# In[3]:


# 函数调用阶段

def load_dat_file(file_path):
    # 打开 dat 文件
    with open(file_path, "rb") as f:
        # 读取文件内容
        data = f.read()

    # 解析字节数据
    result = struct.unpack("<" + "d" * (len(data) // 8), data)

    return result


def save_dat_file(file_path, data):
    # 将数据转换为字节数据
    binary_data = struct.pack("<" + "d" * len(data), *data)

    # 保存为二进制dat文件
    with open(file_path, "wb") as f:
        f.write(binary_data)


def log_transformation(data):
    # 取对数
    log_data = np.log10(data)

    # 计算bin_size
    bin_size = log_data[1] - log_data[0]

    # 计算新的起始点和终点
    new_start = log_data[0] - bin_size / 2
    new_end = log_data[-1] + bin_size / 2

    # 根据新的起始点和终点生成新的对数域数据
    new_data = np.linspace(new_start, new_end, len(energy) + 1)

    new_data = np.power(10, new_data)

    return new_data


# ## en转换代码

# In[5]:


base_dir = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"  # 可以使用windows格式
log_file_path = os.path.join(base_dir, "chose.log")  # 仅仅使用chose文件

# 读取参数文件
with open(log_file_path, "r") as f:
    parameters = f.readlines()[3:]

for parameter in parameters:
    parameter = parameter.strip()
    parameter_folder = parameter
    energy_file_path = os.path.join(base_dir, parameter_folder, "calspec", "en.dat")

    # 检查文件是否存在
    if not os.path.exists(energy_file_path):
        print(f"File not found: {energy_file_path}")
        continue

    energy = load_dat_file(energy_file_path)

    # 对能量数据进行取对数并计算bin_size
    new_energy = log_transformation(energy)

    new_energy_file_path = os.path.join(base_dir, parameter_folder, "calspec", "Xspec_en.dat")

    save_dat_file(new_energy_file_path, new_energy)

    #   print(f"Processed file: {new_energy_file_path}")
    print(len(new_energy))
    print(len(energy))
#   print(type(energy))
#   print(type(new_energy))


# ## flux清理代码

# In[7]:

"""
# 新方法
import os
import numpy as np
import struct
import matplotlib.pyplot as plt

# --- 读取和保存函数 ---
def load_dat_file(file_path):
    with open(file_path, "rb") as f:
        data = f.read()
    num_items = len(data) // 8
    return np.array(struct.unpack("<" + "d" * num_items, data))

def save_dat_file(file_path, data):
    with open(file_path, "wb") as f:
        f.write(struct.pack("<" + "d" * len(data), *data))

# --- 清理孤立非零段：左右 window 个点都为 0 才清理 ---
def clean_isolated_peaks_by_segment(flux, window=3):
    flux = np.array(flux, copy=True)
    n = len(flux)
    mask = flux != 0

    segments = []
    i = 0
    while i < n:
        if mask[i]:
            start = i
            while i + 1 < n and mask[i + 1]:
                i += 1
            end = i
            segments.append((start, end))
        i += 1

    for start, end in segments:
        left_block = flux[max(0, start - window):start]
        right_block = flux[end + 1:min(n, end + 1 + window)]

        if len(left_block) == window and len(right_block) == window:
            if np.all(left_block == 0) and np.all(right_block == 0):
                flux[start:end + 1] = 0

    return flux

# --- 主逻辑 ---
base_dir = "/data/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6_slab/data"
log_file_path = os.path.join(base_dir, "chose.log")

# 指定清理窗口大小
WINDOW = 1

# 读取参数列表（跳过前 3 行）
with open(log_file_path, "r") as f:
    parameters = f.readlines()[3:]

for param in parameters:
    param = param.strip()
    calspec_path = os.path.join(base_dir, param, "calspec")
    flux_path = os.path.join(calspec_path, "flux.dat")
    en_path = os.path.join(calspec_path, "en.dat")
    de_path = os.path.join(calspec_path, "de.dat")

    if not os.path.exists(flux_path) or not os.path.exists(en_path) or not os.path.exists(de_path):
        print(f"跳过缺失文件: {param}")
        continue

    flux = load_dat_file(flux_path)
    en = load_dat_file(en_path)
    de = load_dat_file(de_path)

    flux_clean = clean_isolated_peaks_by_segment(flux, window=WINDOW)

    # 保存结果
    save_path = os.path.join(calspec_path, "flux_clean.dat")
    save_dat_file(save_path, flux_clean)
    print(f"✓ Saved: {save_path}")


# In[27]:

"""
# 旧方法，清洗力度没那么大
import os
import numpy as np
import struct


# --- 读取和保存函数 ---
def load_dat_file(file_path):
    with open(file_path, "rb") as f:
        data = f.read()
    num_items = len(data) // 8
    return np.array(struct.unpack("<" + "d" * num_items, data))


def save_dat_file(file_path, data):
    with open(file_path, "wb") as f:
        f.write(struct.pack("<" + "d" * len(data), *data))


# --- 孤立值清理函数 ---
def clean_isolated_peaks_by_segment(flux, window=3):
    """
    清除被左右全为 0 包围的连续非零片段

    Parameters:
        flux (np.ndarray): 原始 flux 数组
        window (int): 左右窗口大小

    Returns:
        cleaned_flux (np.ndarray): 清洗后的 flux 数组（等长）
    """
    flux = np.array(flux, copy=True)
    n = len(flux)
    mask = flux != 0

    # 找出非零段
    segments = []
    i = 0
    while i < n:
        if mask[i]:
            start = i
            while i + 1 < n and mask[i + 1]:
                i += 1
            end = i
            segments.append((start, end))
        i += 1

    for start, end in segments:
        left_block = flux[max(0, start - window):start]
        right_block = flux[end + 1:min(n, end + 1 + window)]

        # 判断左右各 window 个点是否全为 0
        if len(left_block) == window and len(right_block) == window:
            if np.all(left_block == 0) and np.all(right_block == 0):
                flux[start:end + 1] = 0

    return flux


# --- 主批处理逻辑 ---
base_dir = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"
log_file_path = os.path.join(base_dir, "chose.log")

# 读取参数列表
with open(log_file_path, "r") as f:
    parameters = f.readlines()[3:]

# 遍历每个参数目录
for param in parameters:
    param = param.strip()
    calspec_path = os.path.join(base_dir, param, "calspec")
    flux_path = os.path.join(calspec_path, "flux.dat")

    if not os.path.exists(flux_path):
        print(f"跳过缺失文件: {flux_path}")
        continue

    flux = load_dat_file(flux_path)
    flux_clean = clean_isolated_peaks_by_segment(flux, window=2)

    save_path = os.path.join(calspec_path, "flux_clean.dat")
    save_dat_file(save_path, flux_clean)
    print(f"✓ Saved: {save_path}")

# ## 平滑化代码

# In[9]:


import os
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from skimage.restoration import denoise_tv_chambolle
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d


# === 工具函数 ===
def load_dat_file(file_path):
    """读取二进制 .dat 文件，返回 double 类型的 np.ndarray"""
    with open(file_path, "rb") as f:
        return np.frombuffer(f.read(), dtype="<d")


def save_dat_file(file_path, data):
    """保存 np.ndarray 到 .dat 文件（double 类型）"""
    with open(file_path, "wb") as f:
        f.write(np.array(data, dtype="<d").tobytes())


# === 各种平滑方法 ===
def smooth_tv(flux, weight=0.2):
    return denoise_tv_chambolle(flux, weight=weight)


def smooth_lowess(flux, en, frac=0.06):
    return lowess(flux, en, frac=frac, return_sorted=False)


def smooth_rolling(flux, window=15):
    pad = np.pad(flux, (window // 2, window // 2), mode='edge')
    return np.convolve(pad, np.ones(window) / window, mode='valid')


def smooth_savgol(flux, window=51, polyorder=3):
    return savgol_filter(flux, window_length=window, polyorder=polyorder)


def smooth_gaussian(flux, sigma=4):
    return gaussian_filter1d(flux, sigma=sigma)


# === 批量处理路径 ===
base_dir = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"
log_file = os.path.join(base_dir, "chose.log")

with open(log_file, "r") as f:
    parameters = f.readlines()[3:]  # 跳过前三行空白行

for param in parameters:
    param = param.strip()
    folder = os.path.join(base_dir, param, "calspec")
    if not os.path.isdir(folder):
        print(f"❌ Skipping missing folder: {folder}")
        continue

    try:
        en = load_dat_file(os.path.join(folder, "en.dat"))
        flux = load_dat_file(os.path.join(folder, "flux.dat"))
    except Exception as e:
        print(f"⚠️ Error reading {param}: {e}")
        continue

    # === 各种平滑 ===
    try:
        save_dat_file(os.path.join(folder, "flux_smoothed_lowess.dat"), smooth_lowess(flux, en))
        save_dat_file(os.path.join(folder, "flux_smoothed_rolling.dat"), smooth_rolling(flux))
        save_dat_file(os.path.join(folder, "flux_smoothed_tv.dat"), smooth_tv(flux))
        save_dat_file(os.path.join(folder, "flux_smoothed_savgol.dat"), smooth_savgol(flux))
        save_dat_file(os.path.join(folder, "flux_smoothed_gaussian.dat"), smooth_gaussian(flux))
        print(f"✓ {param}: all 5 smoothings saved.")
    except Exception as e:
        print(f"❌ {param} smoothing error: {e}")

# #### clean之后的平滑

# In[12]:


import os
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from skimage.restoration import denoise_tv_chambolle
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d


# === 工具函数 ===
def load_dat_file(file_path):
    """读取二进制 .dat 文件，返回 double 类型的 np.ndarray"""
    with open(file_path, "rb") as f:
        return np.frombuffer(f.read(), dtype="<d")


def save_dat_file(file_path, data):
    """保存 np.ndarray 到 .dat 文件（double 类型）"""
    with open(file_path, "wb") as f:
        f.write(np.array(data, dtype="<d").tobytes())


# === 各种平滑方法 ===
def smooth_tv(flux, weight=0.2):
    return denoise_tv_chambolle(flux, weight=weight)


def smooth_lowess(flux, en, frac=0.06):
    return lowess(flux, en, frac=frac, return_sorted=False)


def smooth_rolling(flux, window=15):
    pad = np.pad(flux, (window // 2, window // 2), mode='edge')
    return np.convolve(pad, np.ones(window) / window, mode='valid')


def smooth_savgol(flux, window=51, polyorder=3):
    return savgol_filter(flux, window_length=window, polyorder=polyorder)


def smooth_gaussian(flux, sigma=4):
    return gaussian_filter1d(flux, sigma=sigma)


# === 批量处理路径 ===
base_dir = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_7_slab/_exports"
log_file = os.path.join(base_dir, "chose.log")

with open(log_file, "r") as f:
    parameters = f.readlines()[3:]  # 跳过前三行空白行

for param in parameters:
    param = param.strip()
    folder = os.path.join(base_dir, param, "calspec")
    if not os.path.isdir(folder):
        print(f"❌ Skipping missing folder: {folder}")
        continue

    try:
        en = load_dat_file(os.path.join(folder, "en.dat"))
        flux_clean = load_dat_file(os.path.join(folder, "flux_clean.dat"))
    except Exception as e:
        print(f"⚠️ Error reading {param}: {e}")
        continue

    # === 各种平滑 ===
    try:
        save_dat_file(os.path.join(folder, "flux_clean_smoothed_lowess.dat"), smooth_lowess(flux_clean, en))
        save_dat_file(os.path.join(folder, "flux_clean_smoothed_rolling.dat"), smooth_rolling(flux_clean))
        save_dat_file(os.path.join(folder, "flux_clean_smoothed_tv.dat"), smooth_tv(flux_clean))
        save_dat_file(os.path.join(folder, "flux_clean_smoothed_savgol.dat"), smooth_savgol(flux_clean))
        save_dat_file(os.path.join(folder, "flux_clean_smoothed_gaussian.dat"), smooth_gaussian(flux_clean))
        print(f"✓ {param}: flux_clean smoothed results saved.")
    except Exception as e:
        print(f"❌ {param} flux_clean smoothing error: {e}")




