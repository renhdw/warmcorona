
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

def generate_xspec_xcm(output_filename, model_path, parameters, output_dir, additional_commands=None):
    """
    生成一个 XSPEC .xcm 文件

    参数：
    - output_filename: 生成的 xcm 文件名（不包含路径）
    - model_path: 自定义模型的路径（如 warmcom.mod）
    - parameters: 一个列表，包含模型参数值（按你当前模型次序）
    - output_dir: 输出文件的目标目录（可含 ~）
    - additional_commands: 额外的 XSPEC 命令（列表，可选）
    """

    # 展开 ~
    output_dir = os.path.expanduser(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # 写入路径
    output_path = os.path.join(output_dir, output_filename)

    # 参数行：一行一个
    param_lines = "\n".join(str(p) for p in parameters)

    # 额外命令
    extra_cmds = "\n".join(additional_commands) if additional_commands else ""

    # 注意：atable{...} 里放模型的绝对路径/可解析路径
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
""".strip() + "\n"

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(xcm_content)

    print(f"✅ XSPEC xcm 文件已生成: {output_path}")


# ================= 配置区 =================

# 目标目录（Linux 路径；通常放在 ObsID 目录内）
OUTPUT_DIR = "~/data/XMM/Mrk509/0601390401/xcm"

# 两个输出文件名（slab / sphere 各一份）
XCM_SLAB   = "fit_model_Tom_slab_Mrk509.xcm"
XCM_SPHERE = "fit_model_tv_Tom_sphere_Mrk509.xcm"

# 两个模型路径
MODEL_SLAB   = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-0.6_5-25_pure_Tom_slab_107.mod"
MODEL_SPHERE = "/home/hdw/data/monk/plot/warmcorona/model/warmcom/warmcom_0.1-1.0_5-25_smoothed_tv_Tom_new_107.mod"

# 自定义模型参数（与你现有的一致）
params = [
    1,          # zTBabs:nH
    3.469575E-2, # zTBabs:Redshift (frozen)
    3.95E-2,    # TBabs:nH (frozen)
    1.5,        # nthComp:Gamma
    100,        # nthComp:kT_e
    3.00000E-03,# nthComp:kT_bb (frozen)
    0,          # nthComp:inp_type (frozen)
    0.10410067, # nthComp:Redshift (= p2)
    1,          # nthComp:norm
    0.4,        # warmcom:te
    12,         # warmcom:tau
    0.04526328, # warmcom:z (= p2)
    1,          # warmcom:norm
    1.5,        # xillver:gamma (= nthComp:Gamma) p14
    1,          # xillver:Afe (frozen)
    300,        # xillver:Ecut (frozen)
    0,          # xillver:logxi
    0.04526328, # xillver:z (= p2)
    30,         # xillver:Incl
    -1.0,       # xillver:refl_frac (frozen)
    1           # xillver:norm
]

# 额外的 XSPEC 命令（保持你原来的）
extra_cmds = [
    "plot ldata emo re",
    "freeze 3",   # 冻结 TBabs:nH
    "freeze 5",   # 冻结 nthcomp:kTe
    "new 8=2",    # nthComp:Redshift = zTBabs:Redshift
    "new 12=2",   # warmcom:z = zTBabs:Redshift
    "new 14=4",   # xillver:gamma = nthcomp:Gamma
    "new 18=2",   # xillver:z = zTBabs:Redshift  （按你原注释保留）
    "freeze 15",  # 冻结 Afe
    "freeze 16",  # 冻结 Ecut
    "freeze 17",  # 冻结 logxi
    "new 0"
]

if __name__ == "__main__":
    # 生成 slab 版 xcm
    generate_xspec_xcm(
        output_filename=XCM_SLAB,
        model_path=MODEL_SLAB,
        parameters=params,
        output_dir=OUTPUT_DIR,
        additional_commands=extra_cmds
    )

    # 生成 sphere 版 xcm（参数/命令完全一致，仅替换 atable 路径）
    generate_xspec_xcm(
        output_filename=XCM_SPHERE,
        model_path=MODEL_SPHERE,
        parameters=params,
        output_dir=OUTPUT_DIR,
        additional_commands=extra_cmds
    )
