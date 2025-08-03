import numpy as np
import os
import subprocess
import time

# -----------------------------
# ✅ Step 1: 写入原始数据路径到 data_path/Xspec_flux.txt
# -----------------------------
base_dir = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6/data"
data_path_file = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6/data_calspec_1000/Xspec_flux.txt"

te_start = 0.100
te_end = 2.000
te_step = 0.1

tau_start = 10.0
tau_end = 20.0
tau_step = 0.1

te_values = np.arange(te_start, te_end + te_step / 2, te_step)
tau_values = np.arange(tau_start, tau_end + tau_step / 2, tau_step)

lines = []

for te in te_values:
    for tau in tau_values:
        subdir = f"te_{te:.3f}_tau_{tau:.3f}"
        base_path = os.path.join(base_dir, subdir)
        calspec_dir = os.path.join(base_path, "calspec")

        lines.append(f"{base_path}/sphere")
        lines.append(f"{calspec_dir}/flux.dat")
        lines.append(f"{calspec_dir}/en.dat")
        lines.append(f"{calspec_dir}/de.dat")
        lines.append("")  # 每组空一行

# 写入文件
os.makedirs(os.path.dirname(data_path_file), exist_ok=True)
with open(data_path_file, "w") as f:
    f.write("\n".join(lines))

print(f"✅ 原始路径已保存到: {data_path_file}")

# -----------------------------
# ✅ Step 2: 创建输出路径
# -----------------------------
base_dir_calspec = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6/data_calspec_1000"

for te in te_values:
    for tau in tau_values:
        subdir = f"te_{te:.3f}_tau_{tau:.3f}/calspec"
        full_path = os.path.join(base_dir_calspec, subdir)
        os.makedirs(full_path, exist_ok=True)

print("✅ 所有输出路径已创建完成。")


# -----------------------------
# ✅ Step 3: 执行 calspec
# -----------------------------
def run_calspec_program(calspec_folder, sphere_path, calspec_parameter):
    command = f"cd '{calspec_folder}' && /home/hdw/data/monk/monk_for_rhy/bin/calspec '{sphere_path}' {calspec_parameter}"
    result = subprocess.run(["bash", "-c", command], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"❌ calspec error in: {calspec_folder}")
        print(result.stderr)
    else:
        print(f"✅ calspec 执行完成于: {calspec_folder}")
        print(result.stdout)

    time.sleep(0.05)


# ✅ 获取用户输入
calspec_parameters = input("请输入 calspec 的参数（三个数字，用空格分隔）：")
calspec_parameters = calspec_parameters.split()
calspec_parameter = f"{float(calspec_parameters[0])} {float(calspec_parameters[1])} {float(calspec_parameters[2])}"

# ✅ 读取地址文件，每 5 行为一组
with open(data_path_file, "r") as f:
    lines = f.read().strip().split("\n")

group_idx = 0
for i in range(0, len(lines), 5):
    if i + 1 >= len(lines): continue

    sphere_path = lines[i].strip()
    # 以 sphere_path 为基准推导 te_xxx_tau_xxx
    relative = os.path.relpath(sphere_path, base_dir)
    group_name = relative.split("/")[0]
    calspec_folder = os.path.join(base_dir_calspec, group_name, "calspec")

    run_calspec_program(calspec_folder, sphere_path, calspec_parameter)
