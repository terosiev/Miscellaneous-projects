#!/bin/bash
################################################################################
# Script 2: Principal Component Analysis
#
# Performs PCA for population stratification control:
# - Computes principal components on LD-pruned SNPs
# - Generates PCA plots
# - Creates covariate files with PCs for association analysis
#
# Usage: sbatch 02_pca.sh [model_number]
################################################################################

#SBATCH --job-name=gwas_pca
#SBATCH --output=logs/pca/gwas_pca_%A_%a.out
#SBATCH --error=logs/pca/gwas_pca_%A_%a.err
#SBATCH --time=02:00:00
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
echo "Starting PCA for Model ${MODEL}"
echo "==============================================="

# Define directories
QC_OUTPUT="${QC_DIR}/model_${MODEL}"
PCA_OUTPUT="${PCA_DIR}/model_${MODEL}"

# Create output directory
mkdir -p "$PCA_OUTPUT"

# Check input files
if [ ! -f "${QC_OUTPUT}/11_pruned_snps.prune.in" ]; then
    echo "ERROR: Pruned SNP list not found. Run QC first."
    exit 1
fi

## STEP 1: Compute PCA ----
echo "[Step 1] Computing ${N_PCS} principal components..."
plink2 --pfile "$GENOTYPE_FILE" \
       --extract "${QC_OUTPUT}/11_pruned_snps.prune.in" \
       --pca $N_PCS \
       --threads $SLURM_CPUS \
       --out "${PCA_OUTPUT}/pca_results"

## STEP 2: Generate PCA plots ----
echo "[Step 2] Generating PCA plots..."
Rscript "${SCRIPTS_DIR}/r/plot_pca.R" \
        "${PCA_OUTPUT}/pca_results.eigenvec" \
        "${PCA_OUTPUT}/pca_results.eigenval" \
        "${PCA_OUTPUT}/pca_plots"

## STEP 3: Extract PCs as covariates ----
echo "[Step 3] Extracting first ${N_PCS_COVAR} PCs as covariates..."

# Extract FID, IID, and first N_PCS_COVAR PCs
PC_COLS=$(seq 3 $((N_PCS_COVAR + 2)) | tr '\n' ',' | sed 's/,$//')
awk -v cols="1,2,${PC_COLS}" 'BEGIN{split(cols,a,",")} {for(i in a) printf "%s%s", $a[i], (i<length(a)?" ":"\n")}' \
    "${PCA_OUTPUT}/pca_results.eigenvec" > "${PCA_OUTPUT}/pcs_covariates.txt"

## STEP 4: Combine with existing covariates ----
echo "[Step 4] Combining PCs with existing covariates..."

COVAR_FILE="${COVARIATE_PREFIX}_model_${MODEL}.txt"
COMBINED_COVAR="${COVARIATE_DIR}/combined_covariates_model_${MODEL}.txt"

if [ -f "$COVAR_FILE" ]; then
    # Extract PC columns only (without FID/IID)
    PC_ONLY=$(seq 3 $((N_PCS_COVAR + 2)) | tr '\n' ',' | sed 's/,$//')
    awk -v cols="${PC_ONLY}" 'BEGIN{split(cols,a,",")} {for(i in a) printf "%s%s", $a[i], (i<length(a)?" ":"\n")}' \
        "${PCA_OUTPUT}/pca_results.eigenvec" > "${PCA_OUTPUT}/pcs_only.txt"
    
    # Combine side by side
    paste "$COVAR_FILE" "${PCA_OUTPUT}/pcs_only.txt" > "$COMBINED_COVAR"
    
    echo "Combined covariates created: $COMBINED_COVAR"
else
    # No existing covariates, use PCs only
    cp "${PCA_OUTPUT}/pcs_covariates.txt" "$COMBINED_COVAR"
    echo "Warning: No existing covariate file found. Using PCs only."
fi

# Generate summary
echo ""
echo "==============================================="
echo "PCA Summary for Model ${MODEL}"
echo "==============================================="
echo "PCs computed: ${N_PCS}"
echo "PCs used as covariates: ${N_PCS_COVAR}"
echo "Samples: $(tail -n +2 ${PCA_OUTPUT}/pca_results.eigenvec | wc -l)"
echo "Output files:"
echo "  - PCA results: ${PCA_OUTPUT}/pca_results.*"
echo "  - PCA plots: ${PCA_OUTPUT}/pca_plots.pdf"
echo "  - Combined covariates: ${COMBINED_COVAR}"
echo "==============================================="
echo "PCA completed successfully!"
echo "==============================================="
