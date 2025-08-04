import struct
import os
import numpy as np
import matplotlib.pyplot as plt


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
    new_data = np.linspace(new_start, new_end, len(energy)+1)

    new_data = np.power(10, new_data)

    return new_data


base_dir = r"\\wsl.localhost\Ubuntu-22.04\home\hdw\data\monk\plot\control_test\data\warmcom_-50"
log_file_path = os.path.join(base_dir, "chose.log")

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

    print(f"Processed file: {new_energy_file_path}")
    print(len(new_energy))
    print(len(energy))
