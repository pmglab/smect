import os
import pandas as pd
import numpy as np
import anndata as ad
from glob import glob

def convert_spatial_to_h5ad(input_dir, output_subdir="h5ad_output", region_ranges=None):
    """
    将空间转录组文本数据转换为h5ad格式，支持用户指定区域范围
    
    参数:
        input_dir (str): 输入文件目录路径
        output_subdir (str): 输出子目录名称，默认为"h5ad_output"
        region_ranges (list): 区域范围列表，每个元素为字典，包含区域名称和坐标范围
            示例: [{"name": "Layer1", "x_range": (27, 30), "y_range": (27, 30)}]
            如果为None，则使用默认范围
    
    返回:
        str: 输出目录路径
    """
    # 设置默认区域范围
    if region_ranges is None:
        region_ranges = [
            {"name": "Layer1", "x_range": (27, 30), "y_range": (27, 30)},
            {"name": "Layer2", "x_range": None, "y_range": None}  # 默认区域
        ]
    
    # 设置输出路径
    output_dir = os.path.join(input_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有匹配的文件
    file_pattern = os.path.join(input_dir, "gene.spatial.expression.txt.*")
    file_paths = glob(file_pattern)
    
    if not file_paths:
        print(f"在 {input_dir} 中未找到匹配的文件")
        return output_dir
    
    # 处理每个文件
    for file_path in file_paths:
        try:
            adata = process_single_file(file_path, region_ranges)
            save_h5ad_file(adata, file_path, output_dir)
            print(f"成功转换: {os.path.basename(file_path)}")
            
            # 打印区域统计信息
            print_region_statistics(adata, region_ranges)
            
        except Exception as e:
            print(f"处理文件 {file_path} 时出错: {str(e)}")
            continue
    
    print(f"转换完成！结果保存在: {output_dir}")
    return output_dir


def process_single_file(file_path, region_ranges):
    """
    处理单个文件，将其转换为AnnData对象
    
    参数:
        file_path (str): 输入文件路径
        region_ranges (list): 区域范围列表
    
    返回:
        ad.AnnData: 转换后的AnnData对象
    """
    # 读取数据
    df = pd.read_csv(file_path, sep="\t")
    
    # 处理列名和基因名
    df, gene_names = extract_gene_names(df)
    
    # 转置数据并设置基因名为列名
    df = df.T
    df.columns = gene_names
    
    # 解析坐标
    coords = parse_coordinates(df.index)
    
    # 创建观测名（添加C前缀）
    obs_names = create_observation_names(df.index)
    
    # 创建AnnData对象
    adata = create_anndata_object(df, obs_names, gene_names)
    
    # 添加空间坐标
    adata.obsm["spatial"] = coords
    adata.layers["count"] = adata.X
    
    # 添加基于用户指定区域的注释
    add_custom_annotations(adata, coords, region_ranges)
    
    return adata


def add_custom_annotations(adata, coords, region_ranges):
    """
    根据用户指定的区域范围添加注释
    
    参数:
        adata (ad.AnnData): AnnData对象
        coords (np.array): 坐标数组
        region_ranges (list): 区域范围列表
    """
    annotations = []
    
    for coord in coords:
        x, y = coord
        assigned_region = None
        
        # 检查点是否在任何一个指定区域内[1](@ref)
        for region in region_ranges:
            if region["name"] == "Layer2" and region["x_range"] is None:
                # 默认区域，不匹配任何特定区域时使用
                continue
                
            x_range = region["x_range"]
            y_range = region["y_range"]
            
            # 检查点是否在当前区域内[1,5](@ref)
            if (x_range[0] <= x <= x_range[1]) and (y_range[0] <= y <= y_range[1]):
                assigned_region = region["name"]
                break
        
        # 如果没有匹配任何特定区域，则分配到默认区域
        if assigned_region is None:
            # 查找默认区域（通常是最后一个）
            for region in region_ranges:
                if region["x_range"] is None or region["name"] == "Layer2":
                    assigned_region = region["name"]
                    break
        
        annotations.append(assigned_region)
    
    adata.obs["annotation"] = pd.Categorical(annotations)


def is_point_in_region(x, y, region):
    """
    判断点是否在指定区域内[1](@ref)
    
    参数:
        x (float): 点的x坐标
        y (float): 点的y坐标
        region (dict): 区域信息，包含x_range和y_range
    
    返回:
        bool: 点是否在区域内
    """
    if region["x_range"] is None or region["y_range"] is None:
        return False
    
    x_range = region["x_range"]
    y_range = region["y_range"]
    
    return (x_range[0] <= x <= x_range[1]) and (y_range[0] <= y <= y_range[1])


def print_region_statistics(adata, region_ranges):
    """
    打印各区域的统计信息
    
    参数:
        adata (ad.AnnData): AnnData对象
        region_ranges (list): 区域范围列表
    """
    print("区域统计信息:")
    for region in region_ranges:
        region_name = region["name"]
        count = (adata.obs["annotation"] == region_name).sum()
        print(f"  {region_name}: {count} 个点")
        
        # 打印区域范围信息
        if region["x_range"] is not None:
            x_range = region["x_range"]
            y_range = region["y_range"]
            print(f"    范围: x={x_range[0]}-{x_range[1]}, y={y_range[0]}-{y_range[1]}")


# 以下辅助函数保持不变（与之前相同）
def extract_gene_names(df):
    """从数据框中提取基因名称"""
    if df.columns[0] == "gene":
        gene_names = df["gene"].values
        df = df.drop(columns="gene")
    else:
        gene_names = df.iloc[:, 0].values
        df = df.drop(df.columns[0], axis=1)
    
    return df, gene_names


def parse_coordinates(index):
    """从索引字符串解析坐标信息"""
    coords = []
    for idx in index:
        coord_parts = idx.split(":")
        if len(coord_parts) >= 2:
            coords.append([int(coord_parts[0]), int(coord_parts[1])])
        else:
            coords.append([0, 0])
    
    return np.array(coords)


def create_observation_names(index):
    """创建观测名称（添加C前缀）"""
    return ["C" + str(idx) for idx in index]


def create_anndata_object(df, obs_names, gene_names):
    """创建AnnData对象"""
    return ad.AnnData(
        X=df.astype(np.float32).values,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=df.columns)
    )


def save_h5ad_file(adata, file_path, output_dir):
    """保存h5ad文件"""
    base_name = os.path.basename(file_path)
    output_path = os.path.join(output_dir, base_name + ".h5ad")
    adata.write_h5ad(output_path)


# 使用示例
if __name__ == "__main__":
    input_directory = "./Sim_hot30/ExpressionPhenotypeSimulateTask/"
    
    # 示例1：使用默认范围（27<=x<=30和27<=y<=30）
    print("=== 使用默认范围 ===")
    convert_spatial_to_h5ad(input_directory)
    
    # 示例2：使用自定义范围
    print("\n=== 使用自定义范围 ===")
    custom_ranges = [
        {"name": "Hotspot1", "x_range": (25, 28), "y_range": (25, 28)},
        {"name": "Hotspot2", "x_range": (32, 35), "y_range": (32, 35)},
        {"name": "Background", "x_range": None, "y_range": None}  # 默认区域
    ]
    convert_spatial_to_h5ad(input_directory, "h5ad_custom", custom_ranges)
    
    # 示例3：使用多个复杂区域
    print("\n=== 使用多个复杂区域 ===")
    complex_ranges = [
        {"name": "Region_A", "x_range": (20, 25), "y_range": (20, 25)},
        {"name": "Region_B", "x_range": (30, 35), "y_range": (30, 35)},
        {"name": "Region_C", "x_range": (40, 45), "y_range": (40, 45)},
        {"name": "Other", "x_range": None, "y_range": None}
    ]
    convert_spatial_to_h5ad(input_directory, "h5ad_complex", complex_ranges)