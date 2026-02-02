#!/bin/bash
################################################################################
# Script 5: Variant Effect Prediction (VEP)
#
# Annotates lead SNPs with functional information:
# - Gene annotations
# - Consequence predictions
# - SIFT and PolyPhen scores
# - Allele frequencies
# - Canonical transcript information
#
# Requires: VEP installation and cache files
# Usage: sbatch 05_vep.sh [model_number]
################################################################################

#SBATCH --job-name=gwas_vep
#SBATCH --output=logs/vep/gwas_vep_%A_%a.out
#SBATCH --error=logs/vep/gwas_vep_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4

# Source configuration
source scripts/00_config.sh

# Determine model number
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    MODEL=${1:-1}
else
    MODEL=$SLURM_ARRAY_TASK_ID
fi

echo "==============================================="
echo "Starting VEP Annotation for Model ${MODEL}"
echo "==============================================="

# Define directories
CLUMP_OUTPUT="${CLUMP_DIR}/model_${MODEL}"
VEP_OUTPUT="${VEP_DIR}/model_${MODEL}"

# Create output directory
mkdir -p "$VEP_OUTPUT"

# Define input/output files
INPUT_FILE="${CLUMP_OUTPUT}/lead_snps.txt"
OUTPUT_FILE="${VEP_OUTPUT}/vep_annotations.txt"
STATS_FILE="${VEP_OUTPUT}/vep_stats.html"

# Check if input file exists and is not empty
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Lead SNPs file not found: $INPUT_FILE"
    echo "Run clumping step first."
    exit 1
fi

if [ ! -s "$INPUT_FILE" ]; then
    echo "WARNING: No lead SNPs to annotate (file is empty)"
    echo "Creating empty output file..."
    touch "$OUTPUT_FILE"
    exit 0
fi

N_SNPS=$(wc -l < "$INPUT_FILE")
echo "Annotating ${N_SNPS} lead SNPs..."
echo ""

# Check if VEP is available
if ! command -v vep &> /dev/null; then
    echo "ERROR: VEP not found. Please install VEP or load the module."
    echo "Try: module load vep"
    exit 1
fi

## Run VEP annotation ----
echo "[Step 1] Running VEP annotation..."

vep -i "$INPUT_FILE" \
    -o "$OUTPUT_FILE" \
    --species ${VEP_SPECIES} \
    --format id \
    --symbol \
    --cache \
    --dir_cache ${VEP_CACHE_DIR} \
    --force_overwrite \
    --sift b \
    --polyphen b \
    --af \
    --biotype \
    --canonical \
    --tab \
    --stats_file "$STATS_FILE" \
    --verbose

# Check if VEP completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "==============================================="
    echo "VEP Annotation Summary for Model ${MODEL}"
    echo "==============================================="
    echo "Input SNPs: ${N_SNPS}"
    echo "Annotations: $OUTPUT_FILE"
    echo "Statistics: $STATS_FILE"
    echo "==============================================="
    echo "VEP annotation completed successfully!"
    echo "==============================================="
else
    echo "ERROR: VEP annotation failed"
    exit 1
fi

## Extract high-impact variants ----
echo ""
echo "[Step 2] Extracting high-impact variants..."

if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    # Count lines (excluding header and comments)
    N_ANNOTATED=$(grep -v "^#" "$OUTPUT_FILE" | tail -n +2 | wc -l)
    
    # Extract variants with HIGH or MODERATE impact
    grep -v "^#" "$OUTPUT_FILE" | tail -n +2 | \
        awk -F'\t' '$14 == "HIGH" || $14 == "MODERATE"' > "${VEP_OUTPUT}/high_impact_variants.txt"
    
    N_HIGH_IMPACT=$(wc -l < "${VEP_OUTPUT}/high_impact_variants.txt")
    
    echo "Annotated variants: ${N_ANNOTATED}"
    echo "High/moderate impact: ${N_HIGH_IMPACT}"
    
    if [ $N_HIGH_IMPACT -gt 0 ]; then
        echo ""
        echo "High/moderate impact variants saved to:"
        echo "${VEP_OUTPUT}/high_impact_variants.txt"
    fi
fi

echo ""
echo "Annotation complete!"
