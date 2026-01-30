#!/bin/bash
# 完整的数据处理与分析流程脚本
# 包含: 数据预处理 → MAGMA分析 → 基因集生成 → 单细胞数据转换 → scDRS计算

set -e  # 遇到错误立即退出

# 配置变量
BASE_DIR="./Sim_hot30"
SIM_ASSOC_DIR="${BASE_DIR}/simulateS.assoc"
MAGMA_DIR="${BASE_DIR}/magma"
SCDRS_DIR="${BASE_DIR}/scDRS"
RESULT_DIR="${BASE_DIR}/results"
FEATHER_DIR="${BASE_DIR}/ExpressionPhenotypeSimulateTask/gssfeather"

# 创建必要的目录
mkdir -p "${MAGMA_DIR}" "${SCDRS_DIR}" "${RESULT_DIR}"

echo "=== 开始完整数据处理流程 ==="

# ============================================================================
# 步骤1: 数据预处理 - 生成MAGMA输入文件
# ============================================================================
echo "步骤1: 数据预处理..."

if [ ! -d "${SIM_ASSOC_DIR}" ]; then
    echo "错误: 基础目录 ${SIM_ASSOC_DIR} 不存在"
    exit 1
fi

for dir in "${SIM_ASSOC_DIR}"/simulateS.assoc_*_random; do
    if [ -d "$dir" ]; then
        echo "正在处理: $(basename "$dir")"
        
        # 检查输入文件
        if [ ! -f "${dir}/variants.hg38.tsv.gz" ] || [ ! -r "${dir}/variants.hg38.tsv.gz" ]; then
            echo "警告: 无法读取 ${dir}/variants.hg38.tsv.gz，跳过"
            continue
        fi
        
        # 处理文件生成MAGMA格式[1](@ref)
        zcat "${dir}/variants.hg38.tsv.gz" 2>/dev/null | 
        awk '
        BEGIN {
            FS = "\t"
            OFS = "\t"
            print "SNP", "CHR", "BP", "P", "NOBS"
        }
        /^#CHROM/ {
            chrom_col = pos_col = rsid_col = p_col = 0
            for (i = 1; i <= NF; i++) {
                if ($i == "#CHROM") chrom_col = i
                if ($i == "POS") pos_col = i
                if ($i == "ID") rsid_col = i
                if ($i == "disease_Assoc@Logistic_Add_P") p_col = i
            }
            if (chrom_col == 0 || pos_col == 0 || rsid_col == 0) {
                print "错误: 缺少必要的列" > "/dev/stderr"
                exit 1
            }
            next
        }
        !/^#/ && $1 ~ /^[0-9]+$/ {
            snp = (rsid_col && $rsid_col != "") ? $rsid_col : "NA"
            chr = (chrom_col) ? $chrom_col : "NA"
            pos = (pos_col) ? $pos_col : "NA"
            p_value = (p_col && $p_col != "") ? $p_col : "NA"
            
            print snp, chr, pos, p_value, "8000"
        }' > "${dir}/variants.magma.sumstats.tsv"
        
        if [ $? -eq 0 ]; then
            echo "完成: $(basename "$dir")"
        else
            echo "错误: 处理 $(basename "$dir") 时失败"
        fi
    fi
done

echo "数据预处理完成!"

# ============================================================================
# 步骤2: MAGMA分析
# ============================================================================
echo "步骤2: 运行MAGMA分析..."

for dir_path in "${SIM_ASSOC_DIR}"/simulateS.assoc_*_random; do
    dir_name=$(basename "$dir_path")
    number=$(echo "$dir_name" | sed -n 's/.*simulateS.assoc_\([0-9]*\)_random.*/\1/p')
    
    if [ -n "$number" ]; then
        echo "正在处理: $dir_name (编号: $number)"
        
        input_file="$dir_path/variants.magma.sumstats.tsv"
        
        if [ -f "$input_file" ]; then
            # 运行MAGMA命令[1,4](@ref)
            ./magma \
                --bfile ./simLDref/variants.annot.hg38 \
                --pval "$input_file" use='SNP,P' ncol='NOBS' \
                --gene-annot ./simulation/simulation.genes.annot \
                --out "${MAGMA_DIR}/simhot30_s${number}"
            
            if [ $? -eq 0 ]; then
                echo "✓ 成功处理: $dir_name"
            else
                echo "✗ 处理失败: $dir_name"
            fi
        else
            echo "✗ 输入文件不存在: $input_file"
        fi
    else
        echo "✗ 无法从文件夹名提取数字: $dir_name"
    fi
done

echo "MAGMA分析完成!"

# ============================================================================
# 步骤3: 添加基因符号
# ============================================================================
echo "步骤3: 添加基因符号..."

for file in "${MAGMA_DIR}"/*.genes.out; do
    base=$(basename "$file")
    disease=${base%.genes.out}
    out_file="${file%.out}.with_symbol.out"
    
    echo "处理: $base -> ${disease}.genes.with_symbol.out"
    
    awk 'NR==FNR {map[$1]=$2; next} 
         FNR==1 {print "GENE\tSYMBOL\t" $0; next} 
         $1 in map {print $1 "\t" map[$1] "\t" $0}' \
         ./magma/id_to_symbol.txt "$file" > "$out_file"
done

echo "基因符号添加完成!"

# ============================================================================
# 步骤4: 生成Z值汇总表 (Python脚本)
# ============================================================================
echo "步骤4: 生成Z值汇总表..."

python3 << 'EOF'
import pandas as pd
import os
import glob

input_dir = "./Sim_hot30/magma/"
output_dir = "./Sim_hot30/magma/"
os.makedirs(output_dir, exist_ok=True)

file_pattern = os.path.join(input_dir, "*.genes.with_symbol.out")
gene_files = glob.glob(file_pattern)

file_groups = {}
for file_path in gene_files:
    base_name = os.path.basename(file_path)
    s_number = base_name.split("_s")[1].split(".")[0]
    
    if s_number not in file_groups:
        file_groups[s_number] = []
    
    file_groups[s_number].append(file_path)

for s_number, files in file_groups.items():
    df_final = pd.DataFrame()
    
    for file_path in files:
        trait_name = os.path.basename(file_path).split(".")[0]
        
        try:
            df = pd.read_csv(file_path, sep='\s+')
            df_trait = df[['SYMBOL', 'ZSTAT']].rename(columns={'SYMBOL': 'Gene', 'ZSTAT': trait_name})
            
            if df_final.empty:
                df_final = df_trait
            else:
                df_final = df_final.merge(df_trait, on='Gene', how='outer')
        except Exception as e:
            print(f"处理文件 {file_path} 时出错: {e}")
    
    output_file = os.path.join(output_dir, f"all_traits_zscore_s{s_number}.tsv")
    df_final.to_csv(output_file, sep='\t', index=False)
    print(f"成功处理了 {len(files)} 个性状文件，创建汇总文件: {output_file}")

print("Z值汇总表生成完成!")
EOF

# ============================================================================
# 步骤5: 生成基因集文件
# ============================================================================
echo "步骤5: 生成基因集文件..."

for zscore_file in "${MAGMA_DIR}"/all_traits_zscore_s*.tsv; do
    suffix=$(basename "$zscore_file" | sed 's/all_traits_zscore_//' | sed 's/\.tsv//')
    out_file="${MAGMA_DIR}/scdrs_simulation_10kb_${suffix}.gs"
    
    scdrs munge-gs \
        --out-file "$out_file" \
        --zscore-file "$zscore_file" \
        --weight zscore \
        --n-max 1000
    
    echo "处理完成: $zscore_file -> $out_file"
done

echo "基因集文件生成完成!"

# ============================================================================
# 步骤6: 转换单细胞数据格式 (Python脚本)
# ============================================================================
echo "步骤6: 转换单细胞数据格式..."

python3 << 'EOF'
import pandas as pd
import anndata as ad
import os
from glob import glob

feather_dir = "./Sim_hot30/ExpressionPhenotypeSimulateTask/gssfeather/"
output_dir = "./Sim_hot30/scDRS/"

os.makedirs(output_dir, exist_ok=True)
feather_files = glob(os.path.join(feather_dir, "*.feather"))

for feather_path in feather_files:
    df = pd.read_feather(feather_path)
    df.set_index('HUMAN_GENE_SYM', inplace=True)
    data_matrix = df.T
    data_matrix = data_matrix.astype('float32')
    
    adata = ad.AnnData(X=data_matrix.values)
    adata.obs_names = data_matrix.index.tolist()
    adata.var = pd.DataFrame(index=data_matrix.columns.tolist())
    
    base_name = os.path.basename(feather_path).replace('.feather', '.h5ad')
    output_path = os.path.join(output_dir, base_name)
    adata.write_h5ad(output_path)
    print(f"Saved: {output_path}")

print("单细胞数据格式转换完成!")
EOF

# ============================================================================
# 步骤7: 运行scDRS分析
# ============================================================================
echo "步骤7: 运行scDRS分析..."

mkdir -p "${RESULT_DIR}"

for i in {0..99}; do
    H5AD_FILE="${SCDRS_DIR}/gene.spatial.expression.txt.${i}_gene_marker_score.h5ad"
    GS_FILE="${MAGMA_DIR}/scdrs_simulation_10kb_s${i}.gs"
    
    if [ -f "$H5AD_FILE" ]; then
        BASENAME=$(basename "$H5AD_FILE" .h5ad)
        TEMP_OUT_DIR="${RESULT_DIR}/temp_${BASENAME}"
        mkdir -p "$TEMP_OUT_DIR"
        
        echo "正在处理文件: $BASENAME (编号: $i)"
        
        scdrs compute-score \
            --h5ad-file "$H5AD_FILE" \
            --h5ad-species human \
            --gs-file "$GS_FILE" \
            --gs-species human \
            --out-folder "$TEMP_OUT_DIR" \
            --flag-filter-data True \
            --flag-raw-count True \
            --n-ctrl 1000 \
            --flag-return-ctrl-norm-score True
        
        for SCORE_FILE in "$TEMP_OUT_DIR"/*.score.gz; do
            if [ -f "$SCORE_FILE" ]; then
                NEW_NAME=$(basename "$SCORE_FILE" .score.gz)
                mv "$SCORE_FILE" "${RESULT_DIR}/${BASENAME}_${NEW_NAME}.score.gz"
            fi
        done
        
        rm -rf "$TEMP_OUT_DIR"
        echo "已完成处理: $BASENAME"
    else
        echo "文件不存在，跳过: $H5AD_FILE"
    fi
done

echo "=== 完整流程执行完成! ==="
echo "结果保存在: ${RESULT_DIR}"