import pandas as pd
import numpy as np

def reorder_and_sort_para_spec(para_path, spec_path, out_para_path, out_spec_path):
    # 读取原始文件
    para_df = pd.read_csv(para_path)
    spec_df = pd.read_csv(spec_path)

    # 参数列（去掉 No 列）
    param_cols = [c for c in para_df.columns if c != "No"]

    # 先对参数数据做去重排序处理（针对每个参数列独立去重排序）
    # 注意这里我们要保证参数组合的多重排序，而不是单独去重各列（否则会破坏组合对应）
    # 因此这里不做列独立去重，而是保证参数组合唯一且按参数升序排序

    # 去重参数组合（避免重复行）
    para_df_unique = para_df.drop_duplicates(subset=param_cols)

    # 按参数列顺序升序排序整个 DataFrame
    para_df_sorted = para_df_unique.sort_values(by=param_cols).reset_index(drop=True)

    # 将排序好的参数DataFrame索引与原始No对应，重新排序光谱DataFrame
    # 光谱文件按 No 排序，与参数保持一致
    sorted_no_list = para_df_sorted["No"].values

    # 过滤光谱数据，只保留在参数文件中的No
    spec_df_filtered = spec_df[spec_df["No"].isin(sorted_no_list)].copy()

    # 按参数排序后的No顺序对光谱排序，保持对应关系
    spec_df_sorted = spec_df_filtered.set_index("No").loc[sorted_no_list].reset_index()

    # 保存整理后的文件
    para_df_sorted.to_csv(out_para_path, index=False)
    spec_df_sorted.to_csv(out_spec_path, index=False)

    print(f"整理完成！参数文件保存为：{out_para_path}")
    print(f"整理完成！光谱文件保存为：{out_spec_path}")

if __name__ == "__main__":
    base_dir = "~/data/naoc/EOTA"
    para_path = base_dir + "/disk_para.csv"
    spec_path = base_dir + "/disk_spec.csv"
    out_para_path = base_dir + "/disk_para_sorted.csv"
    out_spec_path = base_dir + "/disk_spec_sorted.csv"

    reorder_and_sort_para_spec(para_path, spec_path, out_para_path, out_spec_path)
