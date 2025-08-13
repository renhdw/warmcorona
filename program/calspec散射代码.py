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
        lines.append("")

# 写入文件
os.makedirs(os.path.dirname(data_path_file), exist_ok=True)
with open(data_path_file, "w") as f:
    f.write("\n".join(lines))

print(f"✅ 原始路径已保存到: {data_path_file}")

# -----------------------------
# ✅ Step 2: 创建多个 nsca 输出路径
# -----------------------------
base_dir_calspec = "/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6/data_calspec_1000"

# 获取 nsca 参数列表
multi_nsca = input("是否启用多个 nsca？(y/n): ").strip().lower() == "y"
if multi_nsca:
    nsca_list = input("请输入 nsca 值（例如：0 1 2 -1）: ").strip().split()
    nsca_list = [int(n) for n in nsca_list]
else:
    nsca_list = [None]  # 不指定 nsca，仅执行一次

# 创建输出路径
for te in te_values:
    for tau in tau_values:
        for nsca in nsca_list:
            if nsca is None:
                subdir = f"te_{te:.3f}_tau_{tau:.3f}/calspec_nsca/default"
            else:
                subdir = f"te_{te:.3f}_tau_{tau:.3f}/calspec_nsca/{nsca}"
            full_path = os.path.join(base_dir_calspec, subdir)
            os.makedirs(full_path, exist_ok=True)


print("✅ 所有输出路径已创建完成。")


# -----------------------------
# ✅ Step 3: 执行 calspec（支持 nsca）
# -----------------------------
def run_calspec_program(calspec_folder, sphere_path, base_params, nsca=None):
    """
    运行 calspec 命令，可带 nsca 参数
    """
    nsca_part = f" {nsca}" if nsca is not None else ""
    command = f"cd '{calspec_folder}' && /home/hdw/data/monk/monk_for_rhy/bin/calspec '{sphere_path}' {base_params}{nsca_part}"
    result = subprocess.run(["bash", "-c", command], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"❌ calspec error in: {calspec_folder} (nsca={nsca})")
        print(result.stderr)
    else:
        print(f"✅ calspec 执行完成于: {calspec_folder} (nsca={nsca})")
        print(result.stdout)

    time.sleep(0.01)


# 获取 ne emin emax
base_input = input("请输入 ne emin emax（三个数，用空格分隔）：")
base_params = " ".join([str(float(x)) for x in base_input.split()[:3]])

# 读取地址文件
with open(data_path_file, "r") as f:
    lines = f.read().strip().split("\n")

for i in range(0, len(lines), 5):
    if i + 1 >= len(lines): continue

    sphere_path = lines[i].strip()
    relative = os.path.relpath(sphere_path, base_dir)
    group_name = relative.split("/")[0]

    for nsca in nsca_list:
        if nsca is None:
            calspec_folder = os.path.join(base_dir_calspec, group_name, "calspec_nsca", "default")
        else:
            calspec_folder = os.path.join(base_dir_calspec, group_name, "calspec_nsca", str(nsca))

        run_calspec_program(calspec_folder, sphere_path, base_params, nsca)
