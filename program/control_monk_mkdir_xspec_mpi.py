#!/usr/bin/env python
# coding: utf-8
# 代码用于生成文件夹以及调用monk程序，本版本为任务池版本，可并行运算，仅在linux版本运行。
# In[3]:


import subprocess
import os
import random
from tqdm import tqdm
import time
# from worker import process_point


# In[5]:


def run_wsl_command(command):
    """
    运行 WSL 命令
    """
    subprocess.run(command, shell=True)# 修改命令


# In[7]:


def run_sphere_program(sphere_folder):
    """
    运行 sphere 程序
    """
    command = f"cd '{sphere_folder}' ; sleep 0.01 ; /home/hdw/data/monk/monk_for_rhy/bin/sphere"
    process = subprocess.Popen(command, shell=True) # 修改1命令

    while True:
        if process.poll() is not None:
            print(f"'sphere' program in '{sphere_folder}' has finished.")
            break

    return process.returncode


# In[9]:


def run_calspec_program(calspec_folder, calspec_parameter):
    """
    运行 calspec 程序
    """
    command = f"cd '{calspec_folder}' && /home/hdw/data/monk/monk_for_rhy/bin/calspec ../sphere/ {calspec_parameter}"
    result = subprocess.run([ "bash", "-c", command], capture_output=True, text=True)

    if result.returncode != 0:
        print("calspec error")
        print(result.stderr)
        return

    time.sleep(0.05)

    print(f"'calspec' program in '{calspec_folder}' has finished.")
    print("calspec output:")
    print(result.stdout)


# In[11]:


def choose_points(selection_method, selection_ratio, te_parts, tau_parts, te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound, random_seed=None):
    """
    选择数据点的行或列
    """
    selected_points = []

    if selection_method == 'random':  # 随机选点
        total_points = te_parts * tau_parts
        selected_points_count = round(total_points * selection_ratio)

        if random_seed is not None:
            random.seed(random_seed)

        all_points = [(i, j) for i in range(te_parts) for j in range(tau_parts)]
        selected_points = random.sample(all_points, selected_points_count)

    elif selection_method == 'row':  # 等差选行
        selected_rows = round(te_parts * selection_ratio)
        # 计算最接近选取比例的行数
        selected_rows = int(selected_rows + 0.5) if (selected_rows + 0.5) > selected_rows else int(selected_rows)

        for i in range(selected_rows):
            for j in range(tau_parts):
                selected_points.append((i, j))

    elif selection_method == 'column':  # 等差选列
        selected_columns = round(tau_parts * selection_ratio)
        # 计算最接近选取比例的列数
        selected_columns = int(selected_columns + 0.5) if (selected_columns + 0.5) > selected_columns else int(selected_columns)

        for i in range(te_parts):
            for j in range(selected_columns):
                selected_points.append((i, j))

    elif selection_method == 'all':  # 全选点
        for i in range(te_parts):
            for j in range(tau_parts):
                selected_points.append((i, j))

    # 保存边界点的行和列
    for i in range(te_parts):
        selected_points.append((i, 0))  # 第一列
        selected_points.append((i, tau_parts-1))  # 最后一列
    for j in range(tau_parts):
        selected_points.append((0, j))  # 第一行
        selected_points.append((te_parts-1, j))  # 最后一行

    return selected_points



# In[13]:


# 新加入的代码
def process_point(te_index, tau_index, te_value, tau_value, te_parts, tau_parts,
                  base_file_path_linux, params_path, calspec_parameter,
                  selected_points):

    folder_name = f"te_{te_value:.3f}_tau_{tau_value:.3f}"
    folder = os.path.join(base_file_path_linux , 'data', folder_name).replace('\\', '/') # 注意这里是修改数据存放文件夹的代码
    sphere_folder = os.path.join(folder, 'sphere').replace('\\', '/')
    calspec_folder = os.path.join(folder, 'calspec').replace('\\', '/')
    sphere_params_path = os.path.join(sphere_folder, 'params.txt').replace('\\', '/')

    run_wsl_command(f"mkdir -p '{folder}'")
    run_wsl_command(f"mkdir -p '{sphere_folder}' '{calspec_folder}'")
    run_wsl_command(f"cp '{params_path}' '{sphere_params_path}'")
    run_wsl_command(f"sed -i '2s/.*/te = {te_value:.3f}/' '{sphere_params_path}'")
    run_wsl_command(f"sed -i '3s/.*/tau = {tau_value:.3f}/' '{sphere_params_path}'")

    run_sphere_program(sphere_folder)
    run_calspec_program(calspec_folder, calspec_parameter)

    is_selected = (te_index, tau_index) in selected_points
    if is_selected:
        flux_path = os.path.join(calspec_folder, 'flux.dat').replace('\\', '/')
        en_path = os.path.join(calspec_folder, 'en.dat').replace('\\', '/')
        de_path = os.path.join(calspec_folder, 'de.dat').replace('\\', '/')
        tag = f"te_{te_value:.3f}_tau_{tau_value:.3f}"
        return (flux_path, en_path, de_path, tag)
    return None


# In[15]:


def create_folders(te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound, te_parts, tau_parts,
                   calspec_parameter, selection_method, selection_ratio, random_seed):
    """
    创建文件夹并运行程序
    """
    start_time = time.time()

    te_values = [te_lower_bound + (te_upper_bound - te_lower_bound) * i / (te_parts - 1) for i in range(te_parts)]
    tau_values = [tau_lower_bound + (tau_upper_bound - tau_lower_bound) * i / (tau_parts - 1) for i in range(tau_parts)]

    num_folders = te_parts * tau_parts

    # 选择数据点的行或列
    selected_points = choose_points(selection_method, selection_ratio, te_parts, tau_parts, te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound, random_seed)
    # official_1000_point有10e6光子点。
    # 创建 log 文件

    base_file_path = '/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6/data' # 修改了路径，变为linux格式，本行修改log存储位置
    
    log_file_path = os.path.join(base_file_path , "control.log")
#    log_file_path = r'\\wsl.localhost\Ubuntu-22.04\home\hdw\data\monk\plot\warmcorona\test\test_smooth\data\control.log'
    log_file = open(log_file_path, 'w')

    log_file.write(f"te_lower_bound: {te_lower_bound}\n")
    log_file.write(f"te_upper_bound: {te_upper_bound}\n")
    log_file.write(f"tau_lower_bound: {tau_lower_bound}\n")
    log_file.write(f"tau_upper_bound: {tau_upper_bound}\n")
    log_file.write(f"te_parts: {te_parts}\n")
    log_file.write(f"tau_parts: {tau_parts}\n")
    log_file.write(f"calspec_parameter: {calspec_parameter}\n")
    log_file.write(f"Total_number_of_folders: {num_folders}\n")

    # 添加选点记录的文件
    chose_log_file_path = os.path.join(base_file_path,"chose.log")
#    chose_log_file_path = r'\\wsl.localhost\Ubuntu-22.04\home\hdw\data\monk\plot\warmcorona\test\test_smooth\data\chose.log'
    chose_log_file = open(chose_log_file_path, 'w')

    # 参数位置基础文件
    base_file_path_linux = r'/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6'
    params_path = os.path.join(base_file_path_linux , 'params.txt').replace('\\', '/')

    chose_log_file.write(f"Selection Method: {selection_method}\n")
    chose_log_file.write(f"Selection Ratio: {selection_ratio}\n")
    chose_log_file.write(f"Random Seed: {random_seed}\n")

    # 添加 Xspec_flux 记录的文件
    xspec_flux_log_file_path = os.path.join(base_file_path,"Xspec_flux.log")
#    xspec_flux_log_file_path = r'\\wsl.localhost\Ubuntu-22.04\home\hdw\data\monk\plot\warmcorona\test\test_smooth\data\Xspec_flux.log'
    xspec_flux_log_file = open(xspec_flux_log_file_path, 'w')

    from concurrent.futures import ProcessPoolExecutor, as_completed

    tasks = []
    with ProcessPoolExecutor(max_workers=3) as executor: #用于并行计算
        for i in range(num_folders):
            te_index = i // tau_parts
            tau_index = i % tau_parts
            te_value = te_values[te_index]
            tau_value = tau_values[tau_index]

            future = executor.submit(
                process_point, te_index, tau_index, te_value, tau_value, te_parts, tau_parts,
                base_file_path_linux, params_path, calspec_parameter, selected_points
            )
            tasks.append(future)

        for future in tqdm(as_completed(tasks), total=len(tasks), desc="进度", unit="点"):
            result = future.result()
            if result is not None:
                flux_path, en_path, de_path, tag = result
                xspec_flux_log_file.write(flux_path + '\n')
                xspec_flux_log_file.write(en_path + '\n')
                xspec_flux_log_file.write(de_path + '\n\n')
                chose_log_file.write(tag + '\n')

    log_file.write(f"\n运行时间: {time.time() - start_time}秒\n")
    log_file.close()
    chose_log_file.close()
    xspec_flux_log_file.close()

    print(f"运行时间：{time.time() - start_time} s")


# In[17]:


def main():
    """
    主函数
    """
    test_mode = int(input("请输入是否是测试阶段（1是，0否）："))
    if test_mode == 1:
        te_lower_bound = 0.1
        te_upper_bound = 0.2
        tau_lower_bound = 0.01
        tau_upper_bound = 0.02
        te_parts = 2
        tau_parts = 2
        calspec_parameter = f"{-200} {1e-2} {1e2}"

        # 添加选点功能的示例
        selection_method = input("请选择选点方式（row: 等差选行, column: 等差选列, random: 随机选点, all: 全选点）：")

        if selection_method in ['row', 'column']:  # 等差选点或随机选点
            selection_ratio = float(input("请选择数据点占比（0-1之间）："))
            random_seed = None
            if selection_method == 'random':
                random_seed = int(input("请输入随机种子："))
        elif selection_method == 'all':  # 全选点
            selection_ratio = 1
            random_seed = None

    else:
        te_lower_bound = float(input("请输入te的区间下限："))
        te_upper_bound = float(input("请输入te的区间上限："))
        tau_lower_bound = float(input("请输入tau的区间下限："))
        tau_upper_bound = float(input("请输入tau的区间上限："))
        te_parts = int(input("请输入te的份数："))
        tau_parts = int(input("请输入tau的份数："))
        calspec_parameters = input("请输入calspec的参数（三个数字，用空格分隔）：")
        calspec_parameters = calspec_parameters.split()
        # 将输入的三个数字格式化为字符串
        calspec_parameter = f"{float(calspec_parameters[0])} {float(calspec_parameters[1])} {float(calspec_parameters[2])}"

        selection_method = input("请选择选点方式（row: 等差选行, column: 等差选列, random: 随机选点, all: 全选点）：")

        if selection_method in ['row', 'column']:  # 等差选点或随机选点
            selection_ratio = float(input("请选择数据点占比（0-1之间）："))
            random_seed = None
            if selection_method == 'random':
                random_seed = int(input("请输入随机种子："))
        elif selection_method == 'all':  # 全选点
            selection_ratio = 1
            random_seed = None

    # 创建文件夹并运行程序
    create_folders(te_lower_bound, te_upper_bound, tau_lower_bound, tau_upper_bound, te_parts, tau_parts,
                   calspec_parameter, selection_method, selection_ratio, random_seed)

if __name__ == '__main__':
    main()


# In[ ]:




