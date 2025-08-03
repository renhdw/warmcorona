import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import time

# 读取chose.log文件获取文件夹名称和参数
with open("/home/hdw/data/monk/plot/warmcorona/test/test_smooth/data/chose.log", 'r') as file:
    lines = file.readlines()[3:]  # 跳过前三行
    for line in lines:
        folder_name = line.strip()  # 获取文件夹名称
        te_tau = folder_name.split('_')
        te_value = float(te_tau[1].replace('te', ''))
        tau_value = float(te_tau[3].replace('tau', ''))

        # 构建路径
        base_file_path = f"/home/hdw/data/monk/plot/warmcorona/test/test_smooth"
        data_path = os.path.join(base_file_path, "data", folder_name)
#        data_path = f"/home/hdw/data/monk/plot/warmcorona/test/test_smooth/data/{folder_name}"
        qdp_file = os.path.join(base_file_path, "plot_output" , f"test_{te_value}_{tau_value}.qdp")
        en_file = os.path.join(data_path, "calspec", "en.dat")
        flux_file = os.path.join(data_path, "calspec", "flux.dat")
        de_file = os.path.join(data_path, "calspec", "de.dat")
        qdp_output = os.path.join(data_path, "calspec", f"test_{folder_name}_fluxen_qdp.pdf")
        error_output = os.path.join(data_path, "calspec", f"test_{folder_name}_error_qdp_error.pdf")

        # 创建一个包含xspec命令的脚本文件
        xspec_script_content = f"""
        model atable{{/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom.mod}}
        {te_value}
        {tau_value}
        0
        1
        plot model
        energies 0.01 100 200 log
        setplot command wdata test_{te_value}_{tau_value}.qdp
        setplot command exit
        iplot 
        exit
        """
        print(te_value)
        print(tau_value)
        # 将xspec脚本写入文件
        script_file = f"{folder_name}_xspec_script.xcm"
        with open(script_file, 'w') as file:
            file.write(xspec_script_content)

        # 运行xspec并告诉它从脚本文件执行命令
        subprocess.run(['xspec', '-', script_file])

        # 确保有足够的时间让xspec命令执行
        # time.sleep(0.1)

        # 读取QDP文件
        with open(qdp_file, 'r') as f:
            lines = f.readlines()

        # 提取第一列和第三列的数据
        x_data = []
        y_data = []
        for line in lines:
            if not line.startswith(('!', '@', 'READ')):
                split_line = line.split()
                try:
                    x = float(split_line[0])
                    y = float(split_line[2])
                    x_data.append(x)
                    y_data.append(y)
                except ValueError:
                    continue

        # 加载和处理数据文件
        en = np.fromfile(en_file)
        flux = np.fromfile(flux_file)
        de = np.fromfile(de_file)
        # flux = de * flux
        #  flux /= np.amax(flux)

        # 绘制第一个数据
        fig, ax = plt.subplots()
        ax.loglog(en, flux, label="Monk")

        # 绘制QDP文件的数据
        y_data = np.array(y_data)
       # y_data /= np.amax(y_data)
        ax.loglog(x_data, y_data, label="QDP")
        ax.set_xlabel("Energy [keV]")
        ax.set_ylabel("Normalized Spectrum")
        ax.legend()
        plt.savefig(qdp_output, format='pdf')

        # 计算两个图形的比率
        ratio = np.divide(y_data, flux, out=np.zeros_like(y_data), where=flux!=0)

        # 转换为百分比误差
        percentage_error = (1 - ratio) * 100

        # 创建新的图形来显示误差
        fig, ax = plt.subplots()
        ax.semilogx(en, percentage_error, label="Error")
        ax.set_xlabel("Energy [keV]")
        ax.set_ylabel("Error (%)")
        ax.legend()

        # 在图上标注10个点的数值
        indices = np.linspace(0, len(en) - 1, 10).astype(int)  # 选择10个点
        for i in indices:
            ax.annotate(f'{percentage_error[i]:.2f}', (en[i], percentage_error[i]))
        plt.savefig(error_output, format='pdf')

        # 清除图表以进行下一次迭代
        plt.clf()
