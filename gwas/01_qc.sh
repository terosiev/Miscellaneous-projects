#!/bin/bash
################################################################################
# Script 1: Quality Control
#
# Performs comprehensive quality control on genotype data including:
# - Sample and SNP missingness filtering
# - Sex discrepancy checks
# - Minor allele frequency filtering
# - Hardy-Weinberg equilibrium testing
# - LD pruning for population structure analysis
# - Heterozygosity checks
# - Relatedness filtering
#
# Usage: sbatch 01_qc.sh [model_number]
#        or run with SLURM array: --array=1-N
################################################################################

#SBATCH --job-name=gwas_qc
#SBATCH --output=logs/qc/gwas_qc_%A_%a.out
#SBATCH --error=logs/qc/gwas_qc_%A_%a.err
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
    MODEL=${1:-1}  # Use first argument or default to 1
else
    MODEL=$SLURM_ARRAY_TASK_ID
fi

echo "==============================================="
echo "Starting QC for Model ${MODEL}"
echo "==============================================="

# Create output directory
create_directories $MODEL
OUTPUT_DIR="${QC_DIR}/model_${MODEL}"

# Define input files
INCLUDE_FILE="${INCLUDE_PREFIX}_model_${MODEL}.txt"
SEX_FILE="${SEX_PREFIX}_model_${MODEL}.txt"
PHENO_FILE="${PHENOTYPE_PREFIX}_model_${MODEL}.txt"

# Check input files exist
for file in "$GENOTYPE_FILE.pgen" "$INCLUDE_FILE" "$SEX_FILE" "$PHENO_FILE"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

## STEP 0: Split PAR regions ----
echo "[Step 0] Splitting pseudoautosomal regions..."
plink2 --pfile "$GENOTYPE_FILE" \
       --keep "$INCLUDE_FILE" \
       --split-par b38 \
       --make-pgen \
       --out "$OUTPUT_DIR/00_split_par"

## STEP 1: Update sex and phenotype ----
echo "[Step 1] Updating sex and phenotype information..."
plink2 --pfile "$OUTPUT_DIR/00_split_par" \
       --update-sex "$SEX_FILE" \
       --pheno "$PHENO_FILE" \
       --make-pgen \
       --out "$OUTPUT_DIR/01_updated_sex"

## STEP 2: Check missingness ----
echo "[Step 2] Calculating missingness statistics..."
plink2 --pfile "$OUTPUT_DIR/01_updated_sex" \
       --missing \
       --out "$OUTPUT_DIR/02_missingness"

## STEP 3: Filter SNPs by missingness ----
echo "[Step 3] Filtering SNPs with high missingness (>${GENO_THRESHOLD})..."
plink2 --pfile "$OUTPUT_DIR/01_updated_sex" \
       --geno $GENO_THRESHOLD \
       --make-pgen \
       --out "$OUTPUT_DIR/03_filtered_snps"

## STEP 4: Filter samples by missingness ----
echo "[Step 4] Filtering samples with high missingness (>${MIND_THRESHOLD})..."
plink2 --pfile "$OUTPUT_DIR/03_filtered_snps" \
       --mind $MIND_THRESHOLD \
       --make-pgen \
       --out "$OUTPUT_DIR/04_filtered_samples"

## STEP 5: Sex check ----
echo "[Step 5] Performing sex discrepancy check..."
plink2 --pfile "$OUTPUT_DIR/04_filtered_samples" \
       --check-sex max-female-xf=$MAX_FEMALE_XF min-male-xf=$MIN_MALE_XF \
       --out "$OUTPUT_DIR/05_sex_check"

## STEP 6: Calculate MAF ----
echo "[Step 6] Calculating minor allele frequencies..."
plink2 --pfile "$OUTPUT_DIR/04_filtered_samples" \
       --chr $CHROMOSOMES \
       --freq \
       --geno-counts \
       --out "$OUTPUT_DIR/06_frequencies"

## STEP 7: Filter by MAF ----
echo "[Step 7] Filtering SNPs by MAF (>${MAF_THRESHOLD})..."
plink2 --pfile "$OUTPUT_DIR/04_filtered_samples" \
       --chr $CHROMOSOMES \
       --maf $MAF_THRESHOLD \
       --make-pgen \
       --out "$OUTPUT_DIR/07_maf_filtered"

## STEP 8: Hardy-Weinberg equilibrium test ----
echo "[Step 8] Testing Hardy-Weinberg equilibrium..."
plink2 --pfile "$OUTPUT_DIR/07_maf_filtered" \
       --hardy \
       --out "$OUTPUT_DIR/08_hardy"

## STEP 9: HWE filtering in controls ----
echo "[Step 9] Identifying SNPs in HWE in controls..."
plink2 --pfile "$OUTPUT_DIR/07_maf_filtered" \
       --hwe $HWE_THRESHOLD midp \
       --keep-if ${PHENO_NAME} == 1 \
       --write-snplist \
       --out "$OUTPUT_DIR/09_hwe_pass_controls"

## STEP 10: Apply HWE filter ----
echo "[Step 10] Filtering to SNPs in HWE..."
plink2 --pfile "$OUTPUT_DIR/07_maf_filtered" \
       --extract "$OUTPUT_DIR/09_hwe_pass_controls.snplist" \
       --make-pgen \
       --out "$OUTPUT_DIR/10_hwe_filtered"

## STEP 11: LD pruning for PCA ----
echo "[Step 11] Pruning SNPs for population structure analysis..."
plink2 --pfile "$OUTPUT_DIR/10_hwe_filtered" \
       --indep-pairwise $PRUNE_WINDOW $PRUNE_STEP $PRUNE_R2 \
       --out "$OUTPUT_DIR/11_pruned_snps"

## STEP 12: Heterozygosity check ----
echo "[Step 12] Calculating heterozygosity..."
plink2 --pfile "$OUTPUT_DIR/10_hwe_filtered" \
       --extract "$OUTPUT_DIR/11_pruned_snps.prune.in" \
       --het \
       --out "$OUTPUT_DIR/12_heterozygosity"

## STEP 13: Relatedness check ----
echo "[Step 13] Checking for related individuals..."
plink2 --pfile "$OUTPUT_DIR/10_hwe_filtered" \
       --extract "$OUTPUT_DIR/11_pruned_snps.prune.in" \
       --king-cutoff $KING_CUTOFF \
       --make-king-table \
       --out "$OUTPUT_DIR/13_relatedness"

## STEP 14: Create final QC files ----
echo "[Step 14] Creating final QC dataset..."
plink2 --pfile "$OUTPUT_DIR/10_hwe_filtered" \
       --remove "$OUTPUT_DIR/13_relatedness.king.cutoff.out.id" \
       --make-pgen \
       --out "$OUTPUT_DIR/final_qc"

# Generate summary statistics
echo ""
echo "==============================================="
echo "QC Summary for Model ${MODEL}"
echo "==============================================="
echo "Input samples: $(wc -l < ${INCLUDE_FILE})"
echo "Final samples: $(wc -l < ${OUTPUT_DIR}/final_qc.psam) (after QC)"
echo "SNPs passing QC: $(wc -l < ${OUTPUT_DIR}/final_qc.pvar)"
echo "Pruned SNPs for PCA: $(wc -l < ${OUTPUT_DIR}/11_pruned_snps.prune.in)"
echo "Related individuals removed: $(wc -l < ${OUTPUT_DIR}/13_relatedness.king.cutoff.out.id)"
echo "==============================================="
echo "QC completed successfully!"
echo "Output: ${OUTPUT_DIR}/final_qc"
echo "==============================================="
