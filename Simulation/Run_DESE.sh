#!/bin/bash

BASE_DIR="./"
GSSTXT_DIR="${BASE_DIR}/Sim_hot30/ExpressionPhenotypeSimulateTask/gssfeather"
ASSOC_BASE_DIR="${BASE_DIR}/simulateS.assoc"
OUTPUT_BASE_DIR="${BASE_DIR}/Sim_hot30/DESE"
THREADS=4  # Number of concurrent tasks

# Create temporary pipe for concurrency control
tempfifo="batch_temp.$$"
mkfifo ${tempfifo}
exec 6<>${tempfifo}
rm -f ${tempfifo}

# Create thread locks for semaphore
for ((i=0;i<${THREADS};i++)); do
    echo
done >&6

# Processing function
process_file() {
    local gs_file="$1"
    local file_number=$(echo "${gs_file}" | grep -oE '[0-9]+$')
    
    if [ -z "${file_number}" ]; then
        echo "跳过文件: ${gs_file}"
        return
    fi
    
    local assoc_dir="${ASSOC_BASE_DIR}/simulateS.assoc_${file_number}_random"
    local ref_gty_file="${assoc_dir}/GenerateAnnotationBaseTask/variants.annot.gty.hg38.gtb"
    local sum_file="${assoc_dir}/variants.hg38.tsv.gz"
    local output_prefix="${OUTPUT_BASE_DIR}/simulateS_gss_${file_number}"
    
    if [ ! -f "${ref_gty_file}" ] || [ ! -f "${sum_file}" ]; then
        echo "错误: 输入文件不存在 for ${file_number}"
        return
    fi
    
    echo "开始处理编号: ${file_number}"
    JDK_JAVA_OPTIONS=--add-opens=java.base/java.nio=ALL-UNNAMED \
    java -Xmx48g -jar ./kggsum.jar \
        assoc \
        --ref-gty-file "${ref_gty_file}" \
        refG=hg38 \
        --sum-file "${sum_file}" \
        cp12Cols=CHROM,POS,REF,ALT \
        pbsCols=disease_Assoc@Logistic_Add_P,disease_Assoc@Logistic_Add_Beta,disease_Assoc@Logistic_Add_Beta_SE \
        betaType=1 \
        prevalence=0.01 \
        refG=hg38 \
        --threads 12 \
        --gene-model-database refgene \
        --output "${output_prefix}" \
        --permutation-num 100 \
        --gene-score-file "${gs_file}" \
        ignoreSE=Y \
        --gene-multiple-testing bhfdr \
        --gene-p-cut 0.05
    
    if [ $? -eq 0 ]; then
        echo "完成编号: ${file_number}"
    else
        echo "失败编号: ${file_number}"
    fi
}

# Main loop - process files with numbers from 0 to 99
for file_num in {0..99}; do
    gs_file="${GSSTXT_DIR}/gene.spatial.expression.txt.${file_num}"
    
    # Check if file exists
    if [ ! -f "${gs_file}" ]; then
        echo "文件不存在: ${gs_file}，跳过"
        continue
    fi
    
    read -u6
    {
        process_file "${gs_file}"
        echo >&6
    } &
done

wait
exec 6>&-
echo "所有并行任务完成！"
