#!/bin/bash

sampleNum=100
iterNum=$((sampleNum - 1))

java -Xmx56g -jar ./kgga.jar  \
    simulate \
   --input-gty-file ./sdese.variants.annot.hg38.gtb \
   --xqtl-gene-file ./eqtlgenes.txt \
   --output ./Sim_hot30 \
   --threads 20 \
   --min-obs-rate 0.8 \
   --allele-num 2 \
   --local-maf 0.05~0.5 \
   --r2-cut 0.01 \
   --gene-feature-included 0~14 \
   --xqtl-gene maxQTLProportion=0.4 \
               effectScale=1 \
               heritability=0.2 \
   --phenotype-gene maxDriverNum=250 \
                    heritability=0.15 \
                    prevalence=0.01 \
   --spatial-gene chipSize=100 \
                  hotLocations=30:30 \
                  maxHotExpressionGeneNum=500 \
                  dropoutRate=0.4 \
                  chipNum="$sampleNum" \
   --sampling-group groupNum="$sampleNum" \
   --sampling sampleSize=8000 \
              caseProportion=0.5


 
# Use associative arrays to store raw enrichment values for each condition (array)
declare -A enrichment_values
declare -A significant_count

# Define file path patterns

file_pattern2="./Sim_hot30/ExpressionPhenotypeSimulateTask/gene.spatial.expression.txt."
file_pattern1="./Sim_hot30/ExpressionPhenotypeSimulateTask/phenotype.ped."
 
file1="${file_pattern1}0"
java -Xmx48g -jar ./kgga.jar  \
   annotate \
   --input-gty-file ./sdese.variants.annot.hg38.gtb \
   --ped-file "$file1" \
   --output ./Sim_hot30/simulateS.assoc \
   --threads 32 \
   --seq-ac 1 \
   --hwe 1E-5 \
   --rsid-database dbsnp \
   --output-gty-format PLINK_BED
   
# Iterate over all matching files
INPUT_DIR="./Sim_hot30/ExpressionPhenotypeSimulateTask"
OUTPUT_BASE_DIR="./Sim_hot30/simulateS.assoc"

# Create output directory
mkdir -p "$OUTPUT_BASE_DIR"

# Batch process files from 0 to 99
for i in {0..99}; do
    PED_FILE="$INPUT_DIR/phenotype.ped.$i"
    OUTPUT_DIR="$OUTPUT_BASE_DIR/simulateS.assoc_${i}_random"
    
    # Check if input file exists
    if [ ! -f "$PED_FILE" ]; then
        echo "File $PED_FILE does not exist, skipping"
        continue
    fi
    
    # Create corresponding output directory
    mkdir -p "$OUTPUT_DIR"
    
    echo "Processing: $PED_FILE"
    echo "Output to: $OUTPUT_DIR"
    
    # Execute Java command
    java -Xmx48g -jar ./kgga.jar \
         prune \
       --input-gty-file ./sdese.variants.annot.hg38.gtb \
       --ped-file "$PED_FILE" \
       --output "$OUTPUT_DIR" \
       --threads 20 \
       --seq-ac 1 \
       --hwe 1E-5 \
       --min-obs-rate 0.8 \
       --allele-num 2~4 \
       --local-maf 0.05~0.5 \
       --rsid-database dbsnp \
       --gene-feature-included 0~14 \
       --p-cut 1 \
       --assoc method=logistic-add
    # Check command execution status
    if [ $? -eq 0 ]; then
        echo "Successfully processed: $PED_FILE"
    else
        echo "Failed to process: $PED_FILE"
    fi
    
    echo "----------------------------------------"
done

echo "Batch processing completed!"
