import os

# 用于一键处理数据，仅能在linux下运行
def generate_xspec_xcm(output_filename, model_path, parameters, output_dir, additional_commands=None):
    """
    生成一个 XSPEC .xcm 文件

    参数：
    - output_filename: 生成的 xcm 文件名（不包含路径）
    - model_path: 自定义模型的路径（如 warmcom.mod）
    - parameters: 一个列表，包含模型参数值
    - output_dir: 输出文件的目标目录
    - additional_commands: 额外的 XSPEC 命令（可选）
    """

    # 解析 Linux `~` 为用户 Home 目录
    output_dir = os.path.expanduser(output_dir)

    # 确保目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 生成完整的文件路径
    output_path = os.path.join(output_dir, output_filename)

    # 确保参数是字符串格式
    param_lines = "\n".join(str(p) for p in parameters)

    # 额外命令（如果有的话）
    extra_cmds = "\n".join(additional_commands) if additional_commands else ""

    # XCM 文件内容
    xcm_content = f"""
data PN_spectrum_grp.fits 
cpd /xs
setpl energy
setpl add
ignore 1.8-2.2
ignore **-0.3
ignore 10.0-**
lmod relxill ~/data/monk/plot/warmcorona/model/relxill/relxill_model_v2.3
model zTBabs*TBabs(nthComp + atable{{{model_path}}} + xillver)
{param_lines}
{extra_cmds}
    """
# energies 0.2 30 200 log 预备命令
    # 写入文件
    with open(output_path, "w") as f:
        f.write(xcm_content.strip() + "\n")

    print(f"XSPEC xcm 文件已生成: {output_path}")


# 目标目录（Linux 路径）
output_directory = "~/data/XMM/Mrk509/0601390401/xcm"# 现在一般都在ID文件夹内，先经过SAS数据处理
output_filename = "fit_model_fv_clean_mrk509.xcm"
#output_filename = "~/data/XMM/Mrk509/06013390201/xcm/fit_model_lowess_mrk509.xcm"
model_file = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-2.0_10-20_smoothed_fv_clean.mod"

# 自定义模型参数（按你的需求修改）
params = [
    1,  # zTBabs:nH
    3.469575E-2,  # zTBabs:Redshift (frozen)
    3.95E-2,  # TBabs:nH (frozen)
    1.5,  # nthComp:Gamma
    100,  # nthComp:kT_e
    3.00000E-03,  # nthComp:kT_bb (frozen)
    0,  # nthComp:inp_type (frozen)
    3.469575E-2,  # nthComp:Redshift (= p2)
    1,  # nthComp:norm
    1,  # warmcom:te
    12,  # warmcom:tau
    3.469575E-2,  # warmcom:z (= p2)
    1,  # warmcom:norm
    1.5,  # xillver:gamma (= nthComp:Gamma)p14
    1,  # xillver:Afe (frozen)
    300,  # xillver:Ecut (frozen)
    0,  # xillver:logxi
    3.469575E-2,  # xillver:z (= p2)
    30,  # xillver:Incl
    -1.0,  # xillver:refl_frac (frozen)
    1  # xillver:norm
]

# 额外的 XSPEC 命令（可选）
extra_cmds = [
    "plot ldata emo re",
    "freeze 3",  # 冻结 TBabs:nH
    "freeze 5",  # nthcomp:kte
    "new 8=2",  # nthComp:Redshift = zTBabs:Redshift
    "new 12=2",  # warmcom:z = zTBabs:Redshift
    "new 14=4",  # xillver:gamma = nthcom:gamma
    "new 18=2", # xillver: gamma =
    "freeze 15",  # 冻结Afe
    "freeze 16",  # 冻结 Ecut
    "freeze 17",  # 冻结logxi
    "new 0"
]

# 生成 xcm 文件
generate_xspec_xcm(output_filename, model_file, params, output_directory, extra_cmds)
