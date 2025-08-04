import struct
import os
import numpy as np
from heasp import *

def load_dat_file(file_path):
    # 打开 dat 文件
    with open(file_path, "rb") as f:
        # 读取文件内容
        data = f.read()

    # 解析字节数据
    result = struct.unpack("<" + "d" * (len(data) // 8), data)

    return result

def get_parameters_from_log(log_file_path):
    parameters = []
    with open(log_file_path, "r") as f:
        for line in f.readlines()[3:]:
            param_value = line.strip().split(":")[-1]
            # 解析参数值
            if param_value.startswith("te_") and "_tau_" in param_value:
                te_tau_values = param_value.split("_tau_")
                te_value = float(te_tau_values[0][3:])
                tau_value = float(te_tau_values[1])
                parameters.append(te_value)
                parameters.append(tau_value)
            else:
                print(f"Ignoring invalid parameter value: {param_value}")
    return parameters

test = table()

# set table descriptors and the energy array
test.setModelName("warmcom")
test.setModelUnits("num")
test.setisRedshift(True)
test.setisAdditive(True)
test.setisError(False)

base_dir = r"/home/hdw/data/monk/plot/warmcorona/test/test_smooth_10_6/data_calspec_1000"  #基础数据文件夹

# set up the energies. note that the size is one greater
# than that of the array for the model fluxes
energy_file_path = os.path.join(base_dir, "te_0.100_tau_10.000", "calspec", "Xspec_en.dat")
energy = load_dat_file(energy_file_path)
test.setEnergies(energy)

test.setNumIntParams(2)
test.setNumAddParams(0)

# define first parameter and give it 11 values ranging from
# 0.0 to 2.0 in steps of 0.2.
testpar = tableParameter()
testpar.setName("te")
testpar.setInterpolationMethod(0)# 0是线性，1是对数。
testpar.setInitialValue(0.1)
testpar.setDelta(0.01)
testpar.setMinimum(0.1)
testpar.setBottom(0.1)
testpar.setTop(2)
testpar.setMaximum(2)
testpar.setTabulatedValues(np.linspace(0.1, 2, 20))

# and push it onto the vector of parameters
test.pushParameter(testpar)

# define the second parameter and give it 5 values ranging from
# 4.6 to 5.4 in steps of 0.2.
testpar2 = tableParameter("tau", 0, 10, 0.01, 10, 10, 20, 20)
testpar2.setTabulatedValues(np.linspace(10, 20, 101))

# and push it onto the vector of parameters
test.pushParameter(testpar2)

# now set up the spectra. these are arbitrarily calculated, in a real program
# this step would read a file or call a routine.

log_file_path = os.path.join(base_dir, "chose.log")
parameters = get_parameters_from_log(log_file_path)

# print(parameters)

for i in range(len(parameters) // 2):
    te = parameters[i * 2]
    tau = parameters[i * 2 + 1]
    folder_name = f"te_{te:.3f}_tau_{tau:.3f}"
    file_path = os.path.join(base_dir, folder_name, "calspec", "flux_clean_smoothed_tv.dat")# 只需要flux文件
    de_file_path = os.path.join(base_dir, folder_name, "calspec", "de.dat")
    flux = np.multiply(load_dat_file(file_path), load_dat_file(de_file_path))
    testspec = tableSpectrum()
    testspec.setParameterValues(np.array([te, tau]))
    testspec.setFlux(flux)
    test.pushSpectrum(testspec)

# now write out the table.
tablefile = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-2.0_10-20_smoothed_tv_clean_2.mod"
if os.path.exists(tablefile):
    os.remove(tablefile)
status = test.write(tablefile)
if status != 0:
    print("Failed to write test.mod: status =", status)

print(f"✅ 已保存到: {tablefile}")