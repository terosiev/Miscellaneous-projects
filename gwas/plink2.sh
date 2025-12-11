#!/bin/bash
################################################################################
# GWAS Quality Control and Association Analysis Pipeline
# 
# A complete pipeline for genome-wide association studies including:
# 1. Quality control (QC) filtering
# 2. Population stratification analysis (PCA)
# 3. Association testing with covariates
#
# Designed for SLURM job arrays to process multiple models in parallel
# Note: Claude.ai was used to clean the original script
################################################################################


################################################################################
# SCRIPT 1: Quality Control (QC)
################################################################################
#!/bin/bash
#SBATCH --job-name=gwas_qc
#SBATCH --output=logs/qc/gwas_qc_%A_%a.out
#SBATCH --error=logs/qc/gwas_qc_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4  # Process 4 models in parallel

# Load required modules
module load plink2
module load r/4.4.0

# Define model name based on array task ID
# Each array task processes a different model (e.g., different covariate sets)
MODEL="model_${SLURM_ARRAY_TASK_ID}"

# Define input/output directories
PFILE_DIR="data/raw/plink/genotypes"           # Input genotype files
OUTPUT_DIR="results/qc/${MODEL}"               # QC output directory
INCLUDE_FILE="files/include_${MODEL}.txt"      # Samples to include
SEX_FILE="files/sex_${MODEL}.txt"              # Sex information
PHENO_FILE="files/phenotype_${MODEL}.txt"      # Phenotype data

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Running QC pipeline for $MODEL"
echo "=========================================="


# STEP 0: Split pseudoautosomal regions (PAR)
# The PAR regions on X chromosome need special handling as they behave like autosomes
echo "Step 0: Splitting PAR regions..."
plink2 --pfile "$PFILE_DIR" \
       --keep "$INCLUDE_FILE" \
       --split-par b38 \
       --make-pgen \
       --out "$OUTPUT_DIR/split_par"


# STEP 1: Update sex and phenotype information
# Ensure sex and phenotype data are correctly assigned to samples
echo "Step 1: Updating sex and phenotype..."
plink2 --pfile "$OUTPUT_DIR/split_par" \
       --update-sex "$SEX_FILE" \
       --pheno "$PHENO_FILE" \
       --make-pgen \
       --out "$OUTPUT_DIR/updated_sex"


# STEP 2: Calculate missingness rates
# Generate per-sample and per-SNP missingness statistics
echo "Step 2: Calculating missingness..."
plink2 --pfile "$OUTPUT_DIR/updated_sex" \
       --missing \
       --out "$OUTPUT_DIR/missingness"


# STEP 3: Filter SNPs with high missingness (>10%)
# Remove poor quality SNPs before filtering samples
# This order prevents removal of samples due to poor SNP quality
echo "Step 3: Filtering SNPs with >10% missingness..."
plink2 --pfile "$OUTPUT_DIR/updated_sex" \
       --geno 0.1 \
       --make-pgen \
       --out "$OUTPUT_DIR/filtered_snps"


# STEP 4: Filter samples with high missingness (>10%)
# Remove samples with poor genotyping quality
echo "Step 4: Filtering samples with >10% missingness..."
plink2 --pfile "$OUTPUT_DIR/filtered_snps" \
       --mind 0.1 \
       --make-pgen \
       --out "$OUTPUT_DIR/filtered_samples"


# STEP 5: Sex concordance check
# Identify samples with discordant reported vs genetic sex
# max-female-xf: max X chromosome F coefficient for females
# min-male-xf: min X chromosome F coefficient for males
echo "Step 5: Checking sex concordance..."
plink2 --pfile "$OUTPUT_DIR/filtered_samples" \
       --check-sex max-female-xf=0.3 min-male-xf=0.8 \
       --out "$OUTPUT_DIR/sex_discrepancy"


# STEP 6: Calculate allele frequencies
# Compute minor allele frequencies (MAF) and genotype counts
echo "Step 6: Calculating allele frequencies..."
plink2 --pfile "$OUTPUT_DIR/filtered_samples" \
       --chr 1-23 \
       --freq \
       --geno-counts \
       --out "$OUTPUT_DIR/frequencies"


# STEP 7: Filter SNPs by minor allele frequency
# Remove rare variants (MAF < 1%) to reduce false positives
echo "Step 7: Filtering SNPs with MAF < 0.01..."
plink2 --pfile "$OUTPUT_DIR/filtered_samples" \
       --chr 1-23 \
       --maf 0.01 \
       --make-pgen \
       --out "$OUTPUT_DIR/maf"


# STEP 8: Calculate Hardy-Weinberg equilibrium (HWE) statistics
# HWE deviations can indicate genotyping errors or population stratification
echo "Step 8: Testing Hardy-Weinberg equilibrium..."
plink2 --pfile "$OUTPUT_DIR/maf" \
       --hardy \
       --out "$OUTPUT_DIR/hardy"


# STEP 9: Identify SNPs in HWE in controls
# Test HWE specifically in controls (unaffected individuals)
# Affected samples may violate HWE due to disease-SNP associations
# midp: use mid-p adjustment for better small-sample performance
echo "Step 9: Finding SNPs in HWE in controls (p > 1e-6)..."
plink2 --pfile "$OUTPUT_DIR/maf" \
       --hwe 1e-6 midp \
       --keep-if PHENOTYPE == 1 \
       --write-snplist \
       --out "$OUTPUT_DIR/hwe_pass_controls"


# STEP 10: Keep only SNPs passing HWE test
echo "Step 10: Keeping only HWE-passing SNPs..."
plink2 --pfile "$OUTPUT_DIR/maf" \
       --extract "$OUTPUT_DIR/hwe_pass_controls.snplist" \
       --make-pgen \
       --out "$OUTPUT_DIR/hwe"


# STEP 11: LD pruning for downstream analyses
# Create a pruned SNP set for PCA and relatedness estimation
# Window size: 50 SNPs, shift by 5 SNPs, r² threshold: 0.2
echo "Step 11: LD pruning (r² < 0.2)..."
plink2 --pfile "$OUTPUT_DIR/hwe" \
       --indep-pairwise 50 5 0.2 \
       --out "$OUTPUT_DIR/pruned_snps"


# STEP 12: Check for excess heterozygosity
# Identify samples with unusual heterozygosity rates
# High heterozygosity may indicate sample contamination
# Low heterozygosity may indicate inbreeding
echo "Step 12: Calculating heterozygosity rates..."
plink2 --pfile "$OUTPUT_DIR/hwe" \
       --extract "$OUTPUT_DIR/pruned_snps.prune.in" \
       --het \
       --out "$OUTPUT_DIR/heterozygosity"


# STEP 13: Relatedness check
# Identify and remove related individuals (kinship > 0.2)
# Keep unrelated samples to avoid inflation of test statistics
echo "Step 13: Checking for relatedness (kinship cutoff: 0.2)..."
plink2 --pfile "$OUTPUT_DIR/hwe" \
       --extract "$OUTPUT_DIR/pruned_snps.prune.in" \
       --king-cutoff 0.2 \
       --make-king-table \
       --out "$OUTPUT_DIR/relatedness_check"


# STEP 14: Create final QC-filtered dataset
# Remove related individuals identified in step 13
echo "Step 14: Creating final QC dataset..."
plink2 --pfile "$OUTPUT_DIR/hwe" \
       --remove "$OUTPUT_DIR/relatedness_check.king.cutoff.out.id" \
       --make-pgen \
       --out "$OUTPUT_DIR/final_qc"

echo "=========================================="
echo "QC pipeline completed for $MODEL"
echo "Results: $OUTPUT_DIR/final_qc"
echo "=========================================="


################################################################################
# SCRIPT 2: Population Stratification Analysis (PCA)
################################################################################
#!/bin/bash
#SBATCH --job-name=population_pca
#SBATCH --output=logs/pca/pca_%A_%a.out
#SBATCH --error=logs/pca/pca_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4  # Process 4 models

# Load modules
module load plink2
module load r/4.4.0

# Define model
MODEL="model_${SLURM_ARRAY_TASK_ID}"

# Define directories
PFILE_DIR="data/raw/plink/genotypes"
QC_DIR="results/qc/${MODEL}"
OUTPUT_DIR="results/pca/${MODEL}"
SCRIPTS_DIR="scripts/r"
FILES_DIR="files"
COVARIATE_FILE="$FILES_DIR/covariates_${MODEL}.txt"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Running PCA for $MODEL"
echo "=========================================="


# STEP 1: Compute principal components
# Use LD-pruned SNPs to avoid bias from correlated variants
# Extract 10 PCs to capture population structure
echo "Step 1: Computing PCA (10 components)..."
plink2 --pfile "$PFILE_DIR" \
       --extract "$QC_DIR/pruned_snps.prune.in" \
       --pca 10 \
       --threads 4 \
       --out "$OUTPUT_DIR/pca_results_${MODEL}"


# STEP 2: Visualize PCA results
# Generate scree plot and PC1 vs PC2 scatter plot
echo "Step 2: Plotting PCA results..."
Rscript "$SCRIPTS_DIR/plot_pca.R" \
       "$OUTPUT_DIR/pca_results_${MODEL}.eigenvec" \
       "$OUTPUT_DIR/pca_results_${MODEL}.eigenval" \
       "$OUTPUT_DIR/pca_results_${MODEL}"


# STEP 3: Extract principal components as covariates
# First 7 PCs are commonly used to adjust for population stratification
# Format: IID FID PC1 PC2 PC3 PC4 PC5 PC6 PC7
echo "Step 3: Extracting PC covariates (PCs 1-7)..."
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
    "$OUTPUT_DIR/pca_results_${MODEL}.eigenvec" \
    > "$OUTPUT_DIR/covar_PC_${MODEL}.txt"


# STEP 4: Combine PCs with other covariates
# Merge population structure PCs with clinical/demographic covariates
# This creates the final covariate file for association testing
echo "Step 4: Combining covariates..."
if [[ -f "$COVARIATE_FILE" ]]; then 
    # Use process substitution to merge files column-wise
    # Extracts only PC values (columns 3-9) to avoid duplicate IDs
    paste "$COVARIATE_FILE" \
          <(awk '{print $3, $4, $5, $6, $7, $8, $9}' \
          "$OUTPUT_DIR/pca_results_${MODEL}.eigenvec") \
          > "$FILES_DIR/combined_covariates_${MODEL}.txt"
    echo "Covariates merged successfully."
else
    echo "Warning: Base covariate file not found: $COVARIATE_FILE"
    echo "Using PCs only."
    cp "$OUTPUT_DIR/covar_PC_${MODEL}.txt" \
       "$FILES_DIR/combined_covariates_${MODEL}.txt"
fi

echo "=========================================="
echo "PCA completed for $MODEL"
echo "Results: $OUTPUT_DIR"
echo "=========================================="


################################################################################
# SCRIPT 3: Association Analysis
################################################################################
#!/bin/bash
#SBATCH --job-name=gwas_association
#SBATCH --output=logs/assoc/association_%A_%a.out
#SBATCH --error=logs/assoc/association_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4  # Process 4 models

# Load modules
module load plink2
module load r/4.4.0

# Define model
MODEL="model_${SLURM_ARRAY_TASK_ID}"

# Define directories and files
QC_DIR="results/qc/${MODEL}"
OUTPUT_DIR="results/assoc/${MODEL}"
FILES_DIR="files"
SCRIPTS_DIR="scripts/r"

# Input files
PHENOTYPE_FILE="${FILES_DIR}/phenotype_${MODEL}.txt"
COVARIATE_FILE="${FILES_DIR}/combined_covariates_${MODEL}.txt"
GENOTYPE_FILE="${QC_DIR}/final_qc"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Running association analysis for $MODEL"
echo "=========================================="


# STEP 1: Run association analysis
# Perform logistic regression for case-control GWAS
# --pheno-name: Name of phenotype column to test
# --covar-variance-standardize: Standardize continuous covariates
# --chr 1-22: Autosomes only (exclude X, Y, MT)
# --pfilter 1: Keep all SNPs initially (no p-value filtering)
# --ci 0.95: Calculate 95% confidence intervals for odds ratios
# --glm hide-covar: Don't report covariate effects (only SNP effects)
echo "Step 1: Running GWAS (logistic regression)..."
plink2 --pfile ${GENOTYPE_FILE} \
       --pheno ${PHENOTYPE_FILE} \
       --pheno-name PHENOTYPE \
       --covar ${COVARIATE_FILE} \
       --covar-variance-standardize \
       --chr 1-22 \
       --pfilter 1 \
       --ci 0.95 \
       --threads 4 \
       --glm hide-covar \
       --out "${OUTPUT_DIR}/assoc_results_${MODEL}"


# STEP 2: Generate Manhattan and QQ plots
# Visualize genome-wide association results
# Manhattan plot: shows -log10(p) across chromosomes
# QQ plot: checks for genomic inflation/deflation
echo "Step 2: Creating Manhattan and QQ plots..."
Rscript "${SCRIPTS_DIR}/plot_gwas_results.R" \
       "PHENOTYPE" \
       "$OUTPUT_DIR" \
       "$MODEL"

echo "=========================================="
echo "Association analysis completed for $MODEL"
echo "Results: $OUTPUT_DIR"
echo "=========================================="


################################################################################
# EXAMPLE USAGE
################################################################################

# To run the complete pipeline:
#
# 1. Prepare input files for each model:
#    files/include_model_1.txt          # Sample IDs to include
#    files/sex_model_1.txt              # Sex information (FID IID SEX)
#    files/phenotype_model_1.txt        # Phenotypes (FID IID PHENOTYPE)
#    files/covariates_model_1.txt       # Covariates (FID IID AGE SEX ...)
#
# 2. Submit QC job array:
#    sbatch 01_qc_array.sh
#
# 3. After QC completes, submit PCA job array:
#    sbatch 02_pca_array.sh
#
# 4. After PCA completes, submit association job array:
#    sbatch 03_assoc_array.sh
#
# Each script processes all models in parallel using SLURM job arrays.
# Model-specific parameters can be adjusted in the respective input files.


################################################################################
# QUALITY CONTROL SUMMARY
################################################################################

# This pipeline implements standard GWAS QC procedures:
#
# Sample QC:
# - Remove samples with >10% missing genotypes
# - Check sex concordance (reported vs genetic)
# - Remove samples with excess heterozygosity
# - Remove related individuals (kinship > 0.2)
#
# SNP QC:
# - Remove SNPs with >10% missingness
# - Remove SNPs with MAF < 1%
# - Remove SNPs failing HWE (p < 1e-6 in controls)
# - LD prune for PCA (r² < 0.2)
#
# Association Testing:
# - Adjust for population stratification (7 PCs)
# - Include clinical/demographic covariates
# - Standardize continuous covariates
# - Report 95% confidence intervals
#
# Output files include:
# - QC metrics and filtered genotypes
# - PCA results and plots
# - Association results (per-SNP statistics)
# - Manhattan and QQ plots
