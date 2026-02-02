# GWAS Analysis Pipeline

A comprehensive, modular pipeline for genome-wide association studies (GWAS) using PLINK2 and VEP.

## Overview

This pipeline provides end-to-end GWAS analysis from quality control through variant annotation:

- **Quality Control**: Sample/SNP filtering, relatedness checks, population stratification
- **Association Testing**: Logistic regression with covariate adjustment
- **Post-GWAS**: LD clumping and functional annotation

## Features

✅ Modular design - run individual steps or complete pipeline  
✅ SLURM array job support for multiple models  
✅ Comprehensive QC with detailed logging  
✅ Population stratification control via PCA  
✅ Automated visualization (Manhattan, QQ, PCA plots)  
✅ LD-based clumping for independent signals  
✅ Functional annotation with VEP  
✅ Configurable parameters in single file  

## Directory Structure

```
gwas-pipeline/
├── scripts/
│   ├── 00_config.sh              # Configuration file
│   ├── 01_qc.sh                  # Quality control
│   ├── 02_pca.sh                 # PCA analysis
│   ├── 03_association.sh         # GWAS
│   ├── 04_clumping.sh            # LD clumping
│   ├── 05_vep.sh                 # VEP annotation
│   └── r/                        # R visualization scripts
│       ├── plot_pca.R
│       ├── plot_gwas.R
│       └── calculate_lambda.R
├── data/
│   ├── genotypes/                # PLINK2 genotype files (.pgen/.pvar/.psam)
│   └── covariates/               # Phenotype and covariate files
├── results/
│   ├── qc/                       # QC outputs
│   ├── pca/                      # PCA results
│   ├── association/              # GWAS results
│   ├── clumping/                 # Clumped results
│   └── annotation/               # VEP annotations
├── logs/                         # SLURM log files
├── run_pipeline.sh               # Master runner script
└── README.md                     # This file
```

## Requirements

### Software
- PLINK2 (v2.0+)
- R (v4.0+)
- VEP (Variant Effect Predictor)
- SLURM (for job scheduling - optional)

### R Packages
```r
install.packages(c("tidyverse", "qqman", "data.table"))
```

## Input Data Format

### 1. Genotype Files
PLINK2 binary format (.pgen/.pvar/.psam):
```bash
data/genotypes/genotypes.pgen
data/genotypes/genotypes.pvar
data/genotypes/genotypes.psam
```

### 2. Phenotype File
Tab-separated format with FID, IID, and phenotype columns:
```
FID     IID       PHENOTYPE
FAM001  IND001    2
FAM002  IND002    1
```
- Controls coded as 1
- Cases coded as 2

### 3. Covariate File
Tab-separated with FID, IID, and covariate columns:
```
FID     IID       AGE  SEX
FAM001  IND001    45   1
FAM002  IND002    52   2
```

### 4. Sex File
```
FID     IID       SEX
FAM001  IND001    1
FAM002  IND002    2
```
- Males: 1
- Females: 2

### 5. Include File
List of samples to include in analysis:
```
FID     IID
FAM001  IND001
FAM002  IND002
```

## Configuration

Edit `scripts/00_config.sh` to set paths and parameters:

```bash
# Paths
export PROJECT_DIR="/path/to/your/project"
export GENOTYPE_FILE="${PROJECT_DIR}/data/genotypes/genotypes"

# QC thresholds
export GENO_THRESHOLD=0.1      # SNP missingness (10%)
export MIND_THRESHOLD=0.1      # Sample missingness (10%)
export MAF_THRESHOLD=0.01      # Minor allele frequency (1%)
export HWE_THRESHOLD=1e-6      # Hardy-Weinberg p-value

# Association
export PHENO_NAME="PHENOTYPE"  # Column name in phenotype file

# Clumping
export P1=5e-8                 # Index SNP threshold
export R2=0.1                  # LD threshold
export KB=500                  # Window size (kb)
```

## Usage

### Option 1: Run Complete Pipeline

For a single model:
```bash
./run_pipeline.sh 1
```

For all models:
```bash
./run_pipeline.sh all
```

### Option 2: SLURM Array Jobs

Submit all steps as SLURM array job:
```bash
sbatch scripts/01_qc.sh          # --array=1-4
sbatch scripts/02_pca.sh         # --array=1-4
sbatch scripts/03_association.sh # --array=1-4
sbatch scripts/04_clumping.sh    # --array=1-4
sbatch scripts/05_vep.sh         # --array=1-4
```

### Option 3: Individual Steps

Run steps separately:
```bash
bash scripts/01_qc.sh 1
bash scripts/02_pca.sh 1
bash scripts/03_association.sh 1
bash scripts/04_clumping.sh 1
bash scripts/05_vep.sh 1
```

## Pipeline Steps

### Step 1: Quality Control (01_qc.sh)

Performs comprehensive QC:
- Splits pseudoautosomal regions
- Filters samples/SNPs by missingness
- Sex discrepancy checks
- Minor allele frequency filtering
- Hardy-Weinberg equilibrium testing
- LD pruning for PCA
- Heterozygosity checks
- Relatedness filtering (KING > 0.2)

**Outputs**: `results/qc/model_N/final_qc.*`

### Step 2: PCA (02_pca.sh)

Population stratification control:
- Computes 10 principal components
- Generates scree and scatter plots
- Combines PCs with existing covariates

**Outputs**: 
- PCA results: `results/pca/model_N/pca_results.*`
- Combined covariates: `data/covariates/combined_covariates_model_N.txt`

### Step 3: Association Analysis (03_association.sh)

GWAS using logistic regression:
- Tests each SNP for association with phenotype
- Adjusts for covariates (age, sex, PCs, etc.)
- Generates Manhattan and QQ plots
- Calculates genomic inflation factor (λ)

**Outputs**:
- Results: `results/association/model_N/assoc_results.*`
- Significant SNPs: `results/association/model_N/significant_snps_5e8.txt`
- Plots: `results/association/model_N/manhattan_*.pdf`

### Step 4: LD Clumping (04_clumping.sh)

Identifies independent association signals:
- Clumps SNPs based on LD (r² > 0.1, 500kb window)
- Extracts lead SNPs (P < 5×10⁻⁸)
- Annotates with association statistics

**Outputs**:
- Clumped results: `results/clumping/model_N/clumped_results.clumps`
- Lead SNPs: `results/clumping/model_N/lead_snps.txt`

### Step 5: VEP Annotation (05_vep.sh)

Functional annotation of lead SNPs:
- Gene annotations
- Consequence predictions (missense, synonymous, etc.)
- SIFT and PolyPhen scores
- Allele frequencies
- Canonical transcript information

**Outputs**:
- Annotations: `results/annotation/model_N/vep_annotations.txt`
- High-impact variants: `results/annotation/model_N/high_impact_variants.txt`

## Interpreting Results

### QC Summary
Check `results/qc/model_N/final_qc.log` for:
- Number of samples/SNPs passing QC
- Relatedness exclusions
- Heterozygosity outliers

### Association Results
- **Lambda (λ)**: Should be close to 1.0 (range: 1.0-1.1 acceptable)
- **Manhattan plot**: Look for peaks above genome-wide significance line
- **QQ plot**: Points should follow diagonal; deviation indicates inflation

### Clumping
- Independent loci identified
- Lead SNP for each locus
- Use for downstream fine-mapping or functional follow-up

### VEP Annotations
- Focus on HIGH or MODERATE impact variants
- Check SIFT/PolyPhen predictions for missense variants
- Identify candidate causal genes

## Quality Control Thresholds

### Standard Thresholds (Recommended)
- SNP missingness: < 10%
- Sample missingness: < 10%
- MAF: > 1%
- HWE (controls): P > 1×10⁻⁶
- Relatedness: KING coefficient < 0.2 (2nd degree relatives)

### Adjusting Thresholds
More stringent for smaller samples:
- SNP missingness: < 5%
- MAF: > 5%

More lenient for rare variants:
- MAF: > 0.1%

## Troubleshooting

### Issue: No significant associations

**Possible causes:**
- Underpowered study (too few samples)
- Poor phenotype definition
- Population stratification not adequately controlled

**Solutions:**
- Check sample size and case/control ratio
- Review PCA plots for population structure
- Examine QQ plot for genomic inflation

### Issue: High genomic inflation (λ > 1.1)

**Possible causes:**
- Population stratification
- Cryptic relatedness
- Genotyping errors

**Solutions:**
- Include more PCs as covariates
- Lower relatedness cutoff (KING < 0.15)
- Apply stricter QC filters

### Issue: VEP fails

**Possible causes:**
- VEP cache not installed
- Incorrect cache path

**Solutions:**
- Install VEP cache: `vep_install -a cf -s homo_sapiens -y GRCh38`
- Update `VEP_CACHE_DIR` in config

## Output File Formats

### Association Results (.glm.logistic.hybrid)
```
#CHROM  POS       ID          REF  ALT  A1  TEST  OBS_CT  OR      SE      L95     U95     P
1       752566    rs3094315   G    A    A   ADD   10000   1.05    0.03    0.99    1.11    0.12
```

### Clumping Results (.clumps)
```
CHR  F       SNP         BP        P         TOTAL  NSIG  S05   S01   S001  S0001  SP2
1    1       rs3094315   752566    4.2e-10   15     3     2     5     3     2      rs123,rs456,...
```

### VEP Annotations (.txt)
Tab-separated with columns:
- Uploaded_variation
- Location
- Allele
- Gene
- Feature
- Consequence
- IMPACT
- SYMBOL

## Best Practices

1. **Always run QC first** - Clean data is critical
2. **Check PCA plots** - Look for population outliers
3. **Examine QQ plots** - Assess genomic inflation
4. **Use stringent P-value thresholds** - Genome-wide significance (5×10⁻⁸)
5. **Replicate findings** - Validate in independent cohort
6. **Consider multiple testing** - Account for number of tests

## Citation

If you use this pipeline, please cite:

- **PLINK2**: Chang et al., GigaScience 2015
- **qqman**: Turner, bioRxiv 2014  
- **VEP**: McLaren et al., Genome Biology 2016

## Author

**Tero Sievänen**  

## License

This pipeline is provided for academic and research use.

## Version History

- v1.0 (2025-02): Initial release with full QC, association, and annotation pipeline
