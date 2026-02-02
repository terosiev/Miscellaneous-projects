#!/bin/bash
################################################################################
# Script 3: Association Analysis
#
# Performs genome-wide association analysis using logistic regression:
# - Tests association between SNPs and phenotype
# - Adjusts for covariates (including PCs)
# - Generates Manhattan and QQ plots
#
# Usage: sbatch 03_association.sh [model_number]
################################################################################

#SBATCH --job-name=gwas_assoc
#SBATCH --output=logs/association/gwas_assoc_%A_%a.out
#SBATCH --error=logs/association/gwas_assoc_%A_%a.err
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
echo "Starting Association Analysis for Model ${MODEL}"
echo "==============================================="

# Define directories
QC_OUTPUT="${QC_DIR}/model_${MODEL}"
ASSOC_OUTPUT="${ASSOC_DIR}/model_${MODEL}"

# Create output directory
mkdir -p "$ASSOC_OUTPUT"

# Define input files
GENOTYPE_INPUT="${QC_OUTPUT}/final_qc"
PHENO_FILE="${PHENOTYPE_PREFIX}_model_${MODEL}.txt"
COVAR_FILE="${COVARIATE_DIR}/combined_covariates_model_${MODEL}.txt"

# Check input files exist
for file in "${GENOTYPE_INPUT}.pgen" "$PHENO_FILE" "$COVAR_FILE"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

## STEP 1: Run association analysis ----
echo "[Step 1] Running GWAS with logistic regression..."
plink2 --pfile "${GENOTYPE_INPUT}" \
       --pheno "${PHENO_FILE}" \
       --pheno-name ${PHENO_NAME} \
       --covar "${COVAR_FILE}" \
       --covar-variance-standardize \
       --chr ${CHROMOSOMES} \
       --pfilter ${PFILTER} \
       --ci 0.95 \
       --threads ${SLURM_CPUS} \
       --glm hide-covar \
       --out "${ASSOC_OUTPUT}/assoc_results"

## STEP 2: Generate Manhattan and QQ plots ----
echo "[Step 2] Generating Manhattan and QQ plots..."

# Check if results file exists
RESULT_FILE="${ASSOC_OUTPUT}/assoc_results.${PHENO_NAME}.glm.logistic.hybrid"

if [ ! -f "$RESULT_FILE" ]; then
    echo "ERROR: Association results file not found: $RESULT_FILE"
    exit 1
fi

# Run R plotting script
Rscript "${SCRIPTS_DIR}/r/plot_gwas.R" \
        "$RESULT_FILE" \
        "${ASSOC_OUTPUT}" \
        "model_${MODEL}"

## STEP 3: Extract genome-wide significant SNPs ----
echo "[Step 3] Extracting genome-wide significant SNPs..."

# Extract SNPs with P < 5e-8
awk '$12 < 5e-8 {print $0}' "$RESULT_FILE" > "${ASSOC_OUTPUT}/significant_snps_5e8.txt"

# Extract SNPs with P < 1e-5 (suggestive)
awk '$12 < 1e-5 {print $0}' "$RESULT_FILE" > "${ASSOC_OUTPUT}/suggestive_snps_1e5.txt"

# Count significant SNPs
N_SIG=$(tail -n +2 "${ASSOC_OUTPUT}/significant_snps_5e8.txt" | wc -l)
N_SUG=$(tail -n +2 "${ASSOC_OUTPUT}/suggestive_snps_1e5.txt" | wc -l)

## STEP 4: Calculate genomic inflation factor ----
echo "[Step 4] Calculating genomic inflation factor (lambda)..."
Rscript "${SCRIPTS_DIR}/r/calculate_lambda.R" \
        "$RESULT_FILE" \
        "${ASSOC_OUTPUT}/lambda.txt"

LAMBDA=$(cat "${ASSOC_OUTPUT}/lambda.txt" 2>/dev/null || echo "N/A")

# Generate summary
echo ""
echo "==============================================="
echo "Association Analysis Summary for Model ${MODEL}"
echo "==============================================="
echo "Total SNPs tested: $(tail -n +2 $RESULT_FILE | wc -l)"
echo "Genome-wide significant (P < 5e-8): ${N_SIG}"
echo "Suggestive (P < 1e-5): ${N_SUG}"
echo "Genomic inflation factor (λ): ${LAMBDA}"
echo ""
echo "Output files:"
echo "  - Full results: ${RESULT_FILE}"
echo "  - Significant SNPs: ${ASSOC_OUTPUT}/significant_snps_5e8.txt"
echo "  - Manhattan plot: ${ASSOC_OUTPUT}/manhattan_model_${MODEL}.pdf"
echo "  - QQ plot: ${ASSOC_OUTPUT}/qq_model_${MODEL}.pdf"
echo "==============================================="
echo "Association analysis completed successfully!"
echo "==============================================="
