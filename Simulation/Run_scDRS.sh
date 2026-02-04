#!/bin/bash
# Complete data processing and analysis pipeline
# Includes: Data preprocessing → MAGMA analysis → Gene set generation → Single-cell data conversion → scDRS calculation

set -e  # Exit immediately on error

# Configuration variables
BASE_DIR="./Sim_hot30"
SIM_ASSOC_DIR="${BASE_DIR}/simulateS.assoc"
MAGMA_DIR="${BASE_DIR}/magma"
SCDRS_DIR="${BASE_DIR}/scDRS"
RESULT_DIR="${BASE_DIR}/results"
FEATHER_DIR="${BASE_DIR}/ExpressionPhenotypeSimulateTask/gssfeather"

# Create necessary directories
mkdir -p "${MAGMA_DIR}" "${SCDRS_DIR}" "${RESULT_DIR}"

echo "=== Starting complete data processing pipeline ==="

# ============================================================================
# Step 1: Data preprocessing - Generate MAGMA input files
# ============================================================================
echo "Step 1: Data preprocessing..."

if [ ! -d "${SIM_ASSOC_DIR}" ]; then
    echo "Error: Base directory ${SIM_ASSOC_DIR} does not exist"
    exit 1
fi

for dir in "${SIM_ASSOC_DIR}"/simulateS.assoc_*_random; do
    if [ -d "$dir" ]; then
        echo "Processing: $(basename "$dir")"
        
        # Check input files
        if [ ! -f "${dir}/variants.hg38.tsv.gz" ] || [ ! -r "${dir}/variants.hg38.tsv.gz" ]; then
            echo "Warning: Cannot read ${dir}/variants.hg38.tsv.gz, skipping"
            continue
        fi
        
        # Process files to generate MAGMA format[1](@ref)
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
                print "Error: Missing required columns" > "/dev/stderr"
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
            echo "Completed: $(basename "$dir")"
        else
            echo "Error: Failed to process $(basename "$dir")"
        fi
    fi
done

echo "Data preprocessing completed!"

# ============================================================================
# Step 2: MAGMA analysis
# ============================================================================
echo "Step 2: Running MAGMA analysis..."

for dir_path in "${SIM_ASSOC_DIR}"/simulateS.assoc_*_random; do
    dir_name=$(basename "$dir_path")
    number=$(echo "$dir_name" | sed -n 's/.*simulateS.assoc_\([0-9]*\)_random.*/\1/p')
    
    if [ -n "$number" ]; then
        echo "Processing: $dir_name (Number: $number)"
        
        input_file="$dir_path/variants.magma.sumstats.tsv"
        
        if [ -f "$input_file" ]; then
            # Run MAGMA command[1,5](@ref)
            ./magma \
                --bfile ./simLDref/variants.annot.hg38 \
                --pval "$input_file" use='SNP,P' ncol='NOBS' \
                --gene-annot ./simulation/simulation.genes.annot \
                --out "${MAGMA_DIR}/simhot30_s${number}"
            
            if [ $? -eq 0 ]; then
                echo "✓ Successfully processed: $dir_name"
            else
                echo "✗ Processing failed: $dir_name"
            fi
        else
            echo "✗ Input file does not exist: $input_file"
        fi
    else
        echo "✗ Cannot extract number from folder name: $dir_name"
    fi
done

echo "MAGMA analysis completed!"

# ============================================================================
# Step 3: Add gene symbols
# ============================================================================
echo "Step 3: Adding gene symbols..."

for file in "${MAGMA_DIR}"/*.genes.out; do
    base=$(basename "$file")
    disease=${base%.genes.out}
    out_file="${file%.out}.with_symbol.out"
    
    echo "Processing: $base -> ${disease}.genes.with_symbol.out"
    
    awk 'NR==FNR {map[$1]=$2; next} 
         FNR==1 {print "GENE\tSYMBOL\t" $0; next} 
         $1 in map {print $1 "\t" map[$1] "\t" $0}' \
         ./magma/id_to_symbol.txt "$file" > "$out_file"
done

echo "Gene symbols added successfully!"

# ============================================================================
# Step 4: Generate Z-score summary table (Python script)
# ============================================================================
echo "Step 4: Generating Z-score summary table..."

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
            print(f"Error processing file {file_path}: {e}")
    
    output_file = os.path.join(output_dir, f"all_traits_zscore_s{s_number}.tsv")
    df_final.to_csv(output_file, sep='\t', index=False)
    print(f"Successfully processed {len(files)} trait files, created summary file: {output_file}")

print("Z-score summary table generation completed!")
EOF

# ============================================================================
# Step 5: Generate gene set files
# ============================================================================
echo "Step 5: Generating gene set files..."

for zscore_file in "${MAGMA_DIR}"/all_traits_zscore_s*.tsv; do
    suffix=$(basename "$zscore_file" | sed 's/all_traits_zscore_//' | sed 's/\.tsv//')
    out_file="${MAGMA_DIR}/scdrs_simulation_10kb_${suffix}.gs"
    
    scdrs munge-gs \
        --out-file "$out_file" \
        --zscore-file "$zscore_file" \
        --weight zscore \
        --n-max 1000
    
    echo "Processing completed: $zscore_file -> $out_file"
done

echo "Gene set file generation completed!"

# ============================================================================
# Step 6: Convert single-cell data format (Python script)
# ============================================================================
echo "Step 6: Converting single-cell data format..."

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

print("Single-cell data format conversion completed!")
EOF

# ============================================================================
# Step 7: Run scDRS analysis
# ============================================================================
echo "Step 7: Running scDRS analysis..."

mkdir -p "${RESULT_DIR}"

for i in {0..99}; do
    H5AD_FILE="${SCDRS_DIR}/gene.spatial.expression.txt.${i}_gene_marker_score.h5ad"
    GS_FILE="${MAGMA_DIR}/scdrs_simulation_10kb_s${i}.gs"
    
    if [ -f "$H5AD_FILE" ]; then
        BASENAME=$(basename "$H5AD_FILE" .h5ad)
        TEMP_OUT_DIR="${RESULT_DIR}/temp_${BASENAME}"
        mkdir -p "$TEMP_OUT_DIR"
        
        echo "Processing file: $BASENAME (Number: $i)"
        
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
        echo "Completed processing: $BASENAME"
    else
        echo "File does not exist, skipping: $H5AD_FILE"
    fi
done

echo "=== Complete pipeline execution finished! ==="
echo "Results saved in: ${RESULT_DIR}"
