#!/bin/bash

# Error handling: exit immediately on errors and show line numbers
set -euo pipefail
trap 'echo "Error occurred at line $LINENO, command: $BASH_COMMAND"; exit 1' ERR

# Record start time
echo "Script execution started: $(date)"
START_TIME=$(date +%s)

# Configuration variables
WORKDIR_PARENT="./Sim_hot30/ExpressionPhenotypeSimulateTask/gss"
INPUT_DIR="./Sim_hot30/ExpressionPhenotypeSimulateTask/h5ad_output"
BASE_DIR="./Sim_hot30/simulateS.assoc"

# Create output directory
echo "Creating output directory..."
mkdir -p "$WORKDIR_PARENT"

# Phase 1: Process h5ad files to find latent representations
echo "Starting Phase 1: Processing h5ad files to find latent representations..."
find "$INPUT_DIR" -name "gene.spatial.expression.txt.*.h5ad" | while read H5AD_FILE; do
    FILE_NUM=$(basename "$H5AD_FILE" | grep -oE '[0-9]+' | head -1)
    
    if [ "$FILE_NUM" -ge 0 ] && [ "$FILE_NUM" -le 99 ]; then
        SAMPLE_NAME=$(basename "$H5AD_FILE" .h5ad)
        SAMPLE_WORKDIR="${WORKDIR_PARENT}/${SAMPLE_NAME}"
        mkdir -p "$SAMPLE_WORKDIR"
        
        echo "Processing sample: $SAMPLE_NAME (number: $FILE_NUM)"
        
        gsmap run_find_latent_representations \
            --workdir "$SAMPLE_WORKDIR" \
            --sample_name "$SAMPLE_NAME" \
            --input_hdf5_path "$H5AD_FILE" \
            --annotation 'annotation' \
            --data_layer 'count'
    fi
done

wait
echo "Phase 1 completed: $(date)"

# Phase 2: Latent representation to gene mapping
echo "Starting Phase 2: Latent representation to gene mapping..."
find "$INPUT_DIR" -name "gene.spatial.expression.txt.*.h5ad" | while read H5AD_FILE; do
    FILE_NUM=$(basename "$H5AD_FILE" | grep -oE '[0-9]+' | head -1)
    
    if [ "$FILE_NUM" -ge 0 ] && [ "$FILE_NUM" -le 99 ]; then
        SAMPLE_NAME=$(basename "$H5AD_FILE" .h5ad)
        SAMPLE_WORKDIR="${WORKDIR_PARENT}/${SAMPLE_NAME}"
        
        echo "Processing latent to gene mapping: $SAMPLE_NAME"
        
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
echo "Phase 2 completed: $(date)"

# Phase 3: Generate LDSC format GWAS files
echo "Starting Phase 3: Generating LDSC format GWAS files..."
for dir in ${BASE_DIR}/simulateS.assoc_*_random; do
    if [ -d "$dir" ]; then
        echo "Processing: $(basename $dir)"
        
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
            
            echo "Completed: $(basename $dir)"
        else
            echo "Warning: ${dir}/variants.hg38.tsv.gz does not exist"
        fi
    fi
done

echo "Phase 3 completed: $(date)"

# Phase 4: Generate LD scores (with concurrency control)
echo "Starting Phase 4: Generating LD scores..."
MAX_CONCURRENT=4

# Use named pipe for concurrency control
temp_fifo="/tmp/$$.fifo"
mkfifo $temp_fifo
exec 6<>$temp_fifo
rm -f $temp_fifo

# Initialize tokens
for ((i=1; i<=MAX_CONCURRENT; i++)); do
    echo
done >&6

# Iterate through samples and chromosomes
for SAMPLE_NUM in {0..99}; do
    for CHROM in {1..22}; do
        read -u6
        {
            echo "Starting processing sample $SAMPLE_NUM, chromosome $CHROM ..."
            
            gsmap run_generate_ldscore \
                --workdir "./Sim_hot30/ExpressionPhenotypeSimulateTask/gss/gene.spatial.expression.txt.${SAMPLE_NUM}" \
                --sample_name "gene.spatial.expression.txt.${SAMPLE_NUM}" \
                --chrom $CHROM \
                --bfile_root './simLDref/LDSC_Plink/variants.annot.hg38' \
                --keep_snp_root '/public2/lm_data/ch_snp/hm' \
                --gtf_annotation_file './genome_annotation/gtf/gencode.v46.basic.annotation.gtf' \
                --gene_window_size 10000

            echo "Completed sample $SAMPLE_NUM, chromosome $CHROM"
            echo >&6
        } &
    done
done

wait
exec 6>&-
echo "Phase 4 completed: $(date)"

# Phase 5: Run spatial LDSC analysis (with concurrency control)
echo "Starting Phase 5: Running spatial LDSC analysis..."
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
echo "Phase 5 completed: $(date)"

# Calculate total execution time
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

echo "All tasks completed successfully!"
echo "Results saved in: $WORKDIR_PARENT"
echo "Total execution time: $((ELAPSED_TIME / 60)) minutes $((ELAPSED_TIME % 60)) seconds"
