#!/bin/bash
################################################################################
# Configuration File for GWAS Analysis Pipeline
#
# Define all paths, parameters, and settings here
# Edit this file to match your project structure
################################################################################

## PROJECT STRUCTURE ----

# Base directories
export PROJECT_DIR="/path/to/your/project"
export DATA_DIR="${PROJECT_DIR}/data"
export RESULTS_DIR="${PROJECT_DIR}/results"
export SCRIPTS_DIR="${PROJECT_DIR}/scripts"
export LOGS_DIR="${PROJECT_DIR}/logs"

# Input data directories
export GENOTYPE_DIR="${DATA_DIR}/genotypes"
export COVARIATE_DIR="${DATA_DIR}/covariates"

# Output directories (created automatically)
export QC_DIR="${RESULTS_DIR}/qc"
export PCA_DIR="${RESULTS_DIR}/pca"
export ASSOC_DIR="${RESULTS_DIR}/association"
export CLUMP_DIR="${RESULTS_DIR}/clumping"
export VEP_DIR="${RESULTS_DIR}/annotation"

## INPUT FILES ----

# Genotype files (PLINK2 format)
export GENOTYPE_FILE="${GENOTYPE_DIR}/genotypes"  # Without extension (.pgen/.pvar/.psam)

# Covariate files (per model)
# Format: FID IID covar1 covar2 ...
export PHENOTYPE_PREFIX="${COVARIATE_DIR}/phenotype"      # phenotype_model_1.txt
export COVARIATE_PREFIX="${COVARIATE_DIR}/covariates"     # covariates_model_1.txt
export INCLUDE_PREFIX="${COVARIATE_DIR}/include"          # include_model_1.txt
export SEX_PREFIX="${COVARIATE_DIR}/sex"                  # sex_model_1.txt

## ANALYSIS PARAMETERS ----

# Quality control thresholds
export GENO_THRESHOLD=0.1        # SNP missingness threshold (10%)
export MIND_THRESHOLD=0.1        # Sample missingness threshold (10%)
export MAF_THRESHOLD=0.01        # Minor allele frequency threshold (1%)
export HWE_THRESHOLD=1e-6        # Hardy-Weinberg equilibrium p-value
export KING_CUTOFF=0.2           # Relatedness cutoff (2nd degree)

# Sex check parameters
export MAX_FEMALE_XF=0.3         # Maximum X chromosome F coefficient for females
export MIN_MALE_XF=0.8           # Minimum X chromosome F coefficient for males

# LD pruning parameters
export PRUNE_WINDOW=50           # Window size (kb)
export PRUNE_STEP=5              # Step size (kb)
export PRUNE_R2=0.2              # r² threshold

# PCA parameters
export N_PCS=10                  # Number of principal components to compute
export N_PCS_COVAR=6             # Number of PCs to include as covariates

# Association analysis
export PHENO_NAME="PHENOTYPE"    # Phenotype column name
export PFILTER=1                 # Include all SNPs (no p-value filter)
export CHROMOSOMES="1-22"        # Autosomal chromosomes

# Clumping parameters
export P1=5e-8                   # Index SNP threshold (genome-wide significance)
export P2=0.01                   # Secondary SNP threshold
export R2=0.1                    # LD threshold (r²)
export KB=500                    # Clumping window (kb)

## SLURM PARAMETERS ----

# Resource allocation
export SLURM_TIME="02:00:00"     # Default job time
export SLURM_MEM="16G"           # Default memory
export SLURM_CPUS=4              # Default CPUs
export N_MODELS=4                # Number of models to run

## VEP PARAMETERS ----

# VEP annotation settings
export VEP_SPECIES="homo_sapiens"
export VEP_CACHE_DIR="/path/to/vep/cache"  # Update with your VEP cache location

## MODULE LOADING ----

# Function to load required modules
load_modules() {
    module load plink2 2>/dev/null || echo "Warning: plink2 module not found"
    module load r/4.4.0 2>/dev/null || echo "Warning: R module not found"
}

## UTILITY FUNCTIONS ----

# Create necessary directories
create_directories() {
    local model=$1
    mkdir -p "${QC_DIR}/model_${model}"
    mkdir -p "${PCA_DIR}/model_${model}"
    mkdir -p "${ASSOC_DIR}/model_${model}"
    mkdir -p "${CLUMP_DIR}/model_${model}"
    mkdir -p "${VEP_DIR}/model_${model}"
    mkdir -p "${LOGS_DIR}/qc"
    mkdir -p "${LOGS_DIR}/pca"
    mkdir -p "${LOGS_DIR}/association"
    mkdir -p "${LOGS_DIR}/clumping"
    mkdir -p "${LOGS_DIR}/vep"
}

# Print configuration
print_config() {
    echo "==============================================="
    echo "GWAS Pipeline Configuration"
    echo "==============================================="
    echo "Project directory: ${PROJECT_DIR}"
    echo "Number of models: ${N_MODELS}"
    echo "MAF threshold: ${MAF_THRESHOLD}"
    echo "HWE threshold: ${HWE_THRESHOLD}"
    echo "Clumping P1: ${P1}"
    echo "==============================================="
}

# Export all variables
export -f load_modules
export -f create_directories
export -f print_config
