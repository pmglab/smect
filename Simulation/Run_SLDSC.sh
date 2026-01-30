#!/bin/bash

# 设置错误处理：遇到错误立即退出，并显示行号
set -euo pipefail
trap 'echo "错误发生在行 $LINENO，命令: $BASH_COMMAND"; exit 1' ERR

# 记录开始时间
echo "脚本开始执行: $(date)"
START_TIME=$(date +%s)

# 配置变量
WORKDIR_PARENT="./Sim_hot30/ExpressionPhenotypeSimulateTask/gss"
INPUT_DIR="./Sim_hot30/ExpressionPhenotypeSimulateTask/h5ad_output"
BASE_DIR="./Sim_hot30/simulateS.assoc"

# 创建输出目录
echo "创建输出目录..."
mkdir -p "$WORKDIR_PARENT"

# 第一阶段：处理h5ad文件 - 查找潜在表示
echo "开始第一阶段：处理h5ad文件查找潜在表示..."
find "$INPUT_DIR" -name "gene.spatial.expression.txt.*.h5ad" | while read H5AD_FILE; do
    FILE_NUM=$(basename "$H5AD_FILE" | grep -oE '[0-9]+' | head -1)
    
    if [ "$FILE_NUM" -ge 0 ] && [ "$FILE_NUM" -le 99 ]; then
        SAMPLE_NAME=$(basename "$H5AD_FILE" .h5ad)
        SAMPLE_WORKDIR="${WORKDIR_PARENT}/${SAMPLE_NAME}"
        mkdir -p "$SAMPLE_WORKDIR"
        
        echo "正在处理样本: $SAMPLE_NAME (编号: $FILE_NUM)"
        
        gsmap run_find_latent_representations \
            --workdir "$SAMPLE_WORKDIR" \
            --sample_name "$SAMPLE_NAME" \
            --input_hdf5_path "$H5AD_FILE" \
            --annotation 'annotation' \
            --data_layer 'count'
    fi
done

wait
echo "第一阶段完成: $(date)"

# 第二阶段：潜在表示到基因映射
echo "开始第二阶段：潜在表示到基因映射..."
find "$INPUT_DIR" -name "gene.spatial.expression.txt.*.h5ad" | while read H5AD_FILE; do
    FILE_NUM=$(basename "$H5AD_FILE" | grep -oE '[0-9]+' | head -1)
    
    if [ "$FILE_NUM" -ge 0 ] && [ "$FILE_NUM" -le 99 ]; then
        SAMPLE_NAME=$(basename "$H5AD_FILE" .h5ad)
        SAMPLE_WORKDIR="${WORKDIR_PARENT}/${SAMPLE_NAME}"
        
        echo "正在处理潜在表示到基因映射: $SAMPLE_NAME"
        
        gsmap run_latent_to_gene \
            --workdir "$SAMPLE_WORKDIR" \
            --sample_name "$SAMPLE_NAME" \
            --annotation 'annotation' \
            --latent_representation 'latent_GVAE' \
            --num_neighbour 51 \
            --num_neighbour_spatial 201 
    fi
done

wait
echo "第二阶段完成: $(date)"

# 第三阶段：生成LDSC格式GWAS文件
echo "开始第三阶段：生成LDSC格式GWAS文件..."
for dir in ${BASE_DIR}/simulateS.assoc_*_random; do
    if [ -d "$dir" ]; then
        echo "正在处理: $(basename $dir)"
        
        if [ -f "${dir}/variants.hg38.tsv.gz" ]; then
            zcat "${dir}/variants.hg38.tsv.gz" |
            awk '
            BEGIN {
                OFS = "\t"
                print "SNP", "A1", "A2", "Z", "N"
            }
            /^#CHROM/ {
                for (i=1; i<=NF; i++) {
                    if ($i == "REF") ref_col = i
                    if ($i == "ALT") alt_col = i
                    if ($i == "ID") rsid_col = i
                    if ($i == "disease_Assoc@Logistic_Add_Beta") beta_col = i
                    if ($i == "disease_Assoc@Logistic_Add_Beta_SE") se_col = i
                }
                next
            }
            !/^#/ && $1 ~ /^chr/ {
                a1 = $(alt_col)
                a2 = $(ref_col)
                snp = $(rsid_col)
                beta = $(beta_col)
                se = $(se_col)
                z = (se != 0) ? beta / se : "NA"
                print snp, a1, a2, z, "8000"
            }' |
            gzip > "${dir}/variants.ldsc.anno.sumstats.gz"
            
            echo "完成: $(basename $dir)"
        else
            echo "警告: ${dir}/variants.hg38.tsv.gz 不存在"
        fi
    fi
done

echo "第三阶段完成: $(date)"

# 第四阶段：生成LD分数（并发控制）
echo "开始第四阶段：生成LD分数..."
MAX_CONCURRENT=4

# 使用命名管道控制并发
temp_fifo="/tmp/$$.fifo"
mkfifo $temp_fifo
exec 6<>$temp_fifo
rm -f $temp_fifo

# 初始化令牌
for ((i=1; i<=MAX_CONCURRENT; i++)); do
    echo
done >&6

# 遍历样本和染色体
for SAMPLE_NUM in {0..99}; do
    for CHROM in {1..22}; do
        read -u6
        {
            echo "开始处理样本 $SAMPLE_NUM, 染色体 $CHROM ..."
            
            gsmap run_generate_ldscore \
                --workdir "./Sim_hot30/ExpressionPhenotypeSimulateTask/gss/gene.spatial.expression.txt.${SAMPLE_NUM}" \
                --sample_name "gene.spatial.expression.txt.${SAMPLE_NUM}" \
                --chrom $CHROM \
                --bfile_root './simLDref/LDSC_Plink/variants.annot.hg38' \
                --keep_snp_root '/public2/lm_data/ch_snp/hm' \
                --gtf_annotation_file './genome_annotation/gtf/gencode.v46.basic.annotation.gtf' \
                --gene_window_size 10000

            echo "完成样本 $SAMPLE_NUM, 染色体 $CHROM"
            echo >&6
        } &
    done
done

wait
exec 6>&-
echo "第四阶段完成: $(date)"

# 第五阶段：运行空间LDSC分析（并发控制）
echo "开始第五阶段：运行空间LDSC分析..."
MAX_CONCURRENT_JOBS=8

for SAMPLE_NUM in {0..99}; do
    while [ $(jobs -r | wc -l) -ge $MAX_CONCURRENT_JOBS ]; do
        sleep 5
    done

    gsmap run_spatial_ldsc \
        --workdir "./Sim_hot30/ExpressionPhenotypeSimulateTask/gss/gene.spatial.expression.txt.${SAMPLE_NUM}" \
        --sample_name "gene.spatial.expression.txt.${SAMPLE_NUM}" \
        --trait_name 'sim' \
        --sumstats_file "./Sim_hot30/simulateS.assoc/simulateS.assoc_${SAMPLE_NUM}_random/variants.ldsc.anno.sumstats.gz" \
        --num_processes 16 &
done

wait
echo "第五阶段完成: $(date)"

# 计算总执行时间
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

echo "所有任务处理完成！"
echo "结果保存在: $WORKDIR_PARENT"
echo "总执行时间: $((ELAPSED_TIME / 60)) 分钟 $((ELAPSED_TIME % 60)) 秒"