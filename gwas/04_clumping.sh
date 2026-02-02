#!/bin/bash
################################################################################
# Script 4: LD Clumping
#
# Identifies independent association signals by clumping SNPs in LD:
# - Groups SNPs based on LD structure
# - Identifies lead/index SNPs for each locus
# - Reduces redundancy from correlated SNPs
#
# Usage: sbatch 04_clumping.sh [model_number]
################################################################################

#SBATCH --job-name=gwas_clump
#SBATCH --output=logs/clumping/gwas_clump_%A_%a.out
#SBATCH --error=logs/clumping/gwas_clump_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4

# Source configuration
source scripts/00_config.sh

# Load required modules
load_modules

# Determine model number
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    MODEL=${1:-1}
else
    MODEL=$SLURM_ARRAY_TASK_ID
fi

echo "==============================================="
echo "Starting LD Clumping for Model ${MODEL}"
echo "==============================================="

# Define directories
QC_OUTPUT="${QC_DIR}/model_${MODEL}"
ASSOC_OUTPUT="${ASSOC_DIR}/model_${MODEL}"
CLUMP_OUTPUT="${CLUMP_DIR}/model_${MODEL}"

# Create output directory
mkdir -p "$CLUMP_OUTPUT"

# Define input files
GENOTYPE_INPUT="${QC_OUTPUT}/final_qc"
ASSOC_FILE="${ASSOC_OUTPUT}/assoc_results.${PHENO_NAME}.glm.logistic.hybrid"

# Check input files exist
if [ ! -f "${GENOTYPE_INPUT}.pgen" ]; then
    echo "ERROR: Genotype file not found: ${GENOTYPE_INPUT}.pgen"
    exit 1
fi

if [ ! -f "$ASSOC_FILE" ]; then
    echo "ERROR: Association results not found: $ASSOC_FILE"
    exit 1
fi

# Print clumping parameters
echo ""
echo "Clumping parameters:"
echo "  - Index SNP threshold (P1): ${P1}"
echo "  - Secondary SNP threshold (P2): ${P2}"
echo "  - LD threshold (r²): ${R2}"
echo "  - Window size: ${KB} kb"
echo ""

## STEP 1: Perform clumping ----
echo "[Step 1] Performing LD-based clumping..."
plink2 --pfile "${GENOTYPE_INPUT}" \
       --clump "${ASSOC_FILE}" \
       --clump-field P \
       --clump-p1 ${P1} \
       --clump-p2 ${P2} \
       --clump-r2 ${R2} \
       --clump-kb ${KB} \
       --threads ${SLURM_CPUS} \
       --out "${CLUMP_OUTPUT}/clumped_results"

# Check if clumping produced results
if [ ! -f "${CLUMP_OUTPUT}/clumped_results.clumps" ]; then
    echo "WARNING: No significant clumps found (no SNPs passed P1 threshold)"
    echo "Creating empty output files..."
    touch "${CLUMP_OUTPUT}/lead_snps.txt"
    touch "${CLUMP_OUTPUT}/lead_snps_annotated.txt"
    N_CLUMPS=0
else
    ## STEP 2: Extract lead SNPs ----
    echo "[Step 2] Extracting lead SNPs..."
    
    # Extract just the SNP IDs
    awk 'NR > 1 {print $3}' "${CLUMP_OUTPUT}/clumped_results.clumps" > "${CLUMP_OUTPUT}/lead_snps.txt"
    
    ## STEP 3: Annotate lead SNPs with association stats ----
    echo "[Step 3] Annotating lead SNPs with association statistics..."
    
    # Create header
    head -n 1 "$ASSOC_FILE" > "${CLUMP_OUTPUT}/lead_snps_annotated.txt"
    
    # Extract association results for lead SNPs
    while read snp; do
        grep -w "$snp" "$ASSOC_FILE"
    done < "${CLUMP_OUTPUT}/lead_snps.txt" >> "${CLUMP_OUTPUT}/lead_snps_annotated.txt"
    
    # Count clumps
    N_CLUMPS=$(wc -l < "${CLUMP_OUTPUT}/lead_snps.txt")
fi

## STEP 4: Generate clumping summary ----
echo "[Step 4] Generating clumping summary..."

# Create summary file
cat > "${CLUMP_OUTPUT}/clumping_summary.txt" << EOF
LD Clumping Summary for Model ${MODEL}
========================================

Parameters:
  - Index SNP P-value threshold (P1): ${P1}
  - Secondary SNP P-value threshold (P2): ${P2}
  - LD r² threshold: ${R2}
  - Clumping window: ${KB} kb

Results:
  - Number of independent loci: ${N_CLUMPS}
  - Lead SNPs identified: ${N_CLUMPS}

Output files:
  - Clumping results: ${CLUMP_OUTPUT}/clumped_results.clumps
  - Lead SNP IDs: ${CLUMP_OUTPUT}/lead_snps.txt
  - Annotated lead SNPs: ${CLUMP_OUTPUT}/lead_snps_annotated.txt

EOF

# Print summary to screen
cat "${CLUMP_OUTPUT}/clumping_summary.txt"

echo "==============================================="
echo "LD clumping completed successfully!"
echo "Independent loci identified: ${N_CLUMPS}"
echo "==============================================="
