# GWAS Pipeline Quick Start

Get your GWAS analysis running in 10 minutes!

## 1. Prepare Your Data

Organize your data files:

```bash
gwas-pipeline/
├── data/
│   ├── genotypes/
│   │   ├── genotypes.pgen
│   │   ├── genotypes.pvar
│   │   └── genotypes.psam
│   └── covariates/
│       ├── phenotype_model_1.txt
│       ├── covariates_model_1.txt
│       ├── sex_model_1.txt
│       └── include_model_1.txt
```

### Required File Formats

**Phenotype file** (phenotype_model_1.txt):
```
FID     IID       PHENOTYPE
FAM001  IND001    2
FAM002  IND002    1
```

**Covariate file** (covariates_model_1.txt):
```
FID     IID       AGE  SEX
FAM001  IND001    45   1
FAM002  IND002    52   2
```

**Sex file** (sex_model_1.txt):
```
FID     IID       SEX
FAM001  IND001    1
FAM002  IND002    2
```

**Include file** (include_model_1.txt):
```
FID     IID
FAM001  IND001
FAM002  IND002
```

## 2. Configure Pipeline

Edit `scripts/00_config.sh`:

```bash
# Update these paths
export PROJECT_DIR="/your/project/path"
export GENOTYPE_FILE="${PROJECT_DIR}/data/genotypes/genotypes"
export PHENO_NAME="PHENOTYPE"  # Column name in phenotype file
```

## 3. Run Analysis

### Option A: Run Everything
```bash
./run_pipeline.sh 1
```

### Option B: Submit SLURM Jobs
```bash
sbatch scripts/01_qc.sh
sbatch scripts/02_pca.sh
sbatch scripts/03_association.sh
sbatch scripts/04_clumping.sh
sbatch scripts/05_vep.sh
```

## 4. Check Results

```bash
# QC summary
cat results/qc/model_1/final_qc.log

# Association results
head results/association/model_1/significant_snps_5e8.txt

# View Manhattan plot
open results/association/model_1/manhattan_model_1.pdf

# Lead SNPs
head results/clumping/model_1/lead_snps.txt
```

## Expected Runtime

- QC: ~30-60 min (depends on dataset size)
- PCA: ~5-10 min
- Association: ~30-120 min (depends on # SNPs)
- Clumping: ~10-20 min
- VEP: ~5-15 min

**Total: 1-3 hours for typical GWAS**

## Key Output Files

```
results/
├── qc/model_1/final_qc.*              # QC'd genotypes
├── pca/model_1/pca_plots_scatter.pdf  # Population structure
├── association/model_1/
│   ├── assoc_results.*.glm.logistic.hybrid  # Full results
│   ├── manhattan_model_1.pdf          # Manhattan plot
│   ├── qq_model_1.pdf                 # QQ plot
│   └── significant_snps_5e8.txt       # Genome-wide hits
├── clumping/model_1/
│   ├── lead_snps.txt                  # Independent signals
│   └── lead_snps_annotated.txt        # With P-values
└── annotation/model_1/
    └── vep_annotations.txt            # Functional info
```

## Troubleshooting

### "File not found" errors
- Check paths in `00_config.sh`
- Verify file names match expected pattern

### No significant associations
- This is normal for small samples or weak effects
- Check QC metrics and sample size

### VEP fails
- Install VEP cache if not available
- Update `VEP_CACHE_DIR` in config

## Next Steps

1. Review QC plots and metrics
2. Examine Manhattan and QQ plots
3. Look up lead SNPs in GWAS Catalog
4. Plan replication study
5. Consider fine-mapping significant loci

## Need Help?

See full [README.md](README.md) for detailed documentation.
