# Quick Start Guide

Get started with the RNA-seq pipeline in 5 minutes!

## 1. Install Dependencies

```r
# Install required packages (run once)
source("scripts/install_packages.R")
```

## 2. Prepare Your Data

Place your files in the `data/` directory:

```
data/
├── gene_counts.tsv        # Count matrix
└── sample_metadata.tsv    # Sample information
```

### Count Matrix Format

```
gene_id             gene_name    Sample1  Sample2  Sample3
ENSMUSG00000000001  Gnai3        500      450      520
ENSMUSG00000000003  Pbsn         120      130      110
```

### Metadata Format

```
sample     condition
Sample1    Control
Sample2    Control
Sample3    Treatment
```

## 3. Configure Analysis

Edit `scripts/00_config.R` if needed:

```r
# For human data, change:
ORGDB <- "org.Hs.eg.db"
ORGANISM <- "human"
KEGG_ORGANISM <- "hsa"
```

## 4. Run Pipeline

```bash
Rscript run_pipeline.R
```

## 5. Check Results

Results are organized in `results/`:

```
results/
├── qc/     # Quality control plots
├── de/     # Differential expression results
├── ora/    # Pathway enrichment (ORA)
└── gsea/   # Gene set enrichment (GSEA)
```

## Key Output Files

- `results/de/de_summary.csv` - Summary statistics
- `results/de/*_sig_genes.csv` - Significant genes
- `results/de/*_volcano.pdf` - Volcano plots
- `results/de/*_heatmap_top50.pdf` - Heatmaps

## Next Steps

1. Review QC plots in `results/qc/`
2. Check DE summary table
3. Explore enrichment results
4. Customize contrasts or parameters as needed

## Need Help?

See the full [README.md](README.md) for detailed documentation.
