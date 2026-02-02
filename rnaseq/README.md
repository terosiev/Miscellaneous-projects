# RNA-seq Differential Expression Analysis Pipeline

A modular, reproducible pipeline for RNA-seq differential expression analysis using DESeq2, with comprehensive quality control and functional enrichment analysis.

## Overview

This pipeline provides end-to-end analysis of bulk RNA-seq data, including:

- **Quality Control**: Library size assessment, count distribution, PCA, sample correlation
- **Differential Expression**: DESeq2-based analysis with log fold change shrinkage
- **Functional Enrichment**: Over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA) using GO, Reactome, and KEGG databases

## Features

- ✅ Modular design with separate scripts for each analysis step
- ✅ Centralized configuration for easy parameter adjustment
- ✅ Comprehensive QC visualizations (PCA, heatmaps, dispersion plots)
- ✅ Multiple pathway databases (GO, Reactome, KEGG)
- ✅ Both ORA and GSEA approaches for enrichment analysis
- ✅ Automated volcano plots, MA plots, and heatmaps
- ✅ Support for multiple organisms (mouse, human, fly, etc.)

## Directory Structure

```
rnaseq-pipeline/
├── data/                           # Input data directory
│   ├── gene_counts.tsv            # Raw count matrix (required)
│   └── sample_metadata.tsv        # Sample information (required)
├── scripts/                        # Analysis scripts
│   ├── 00_config.R                # Configuration file
│   ├── 01_qc_filtering.R          # Quality control and filtering
│   ├── 02_differential_expression.R  # DE analysis
│   └── 03_functional_enrichment.R    # Pathway enrichment
├── results/                        # Output directory (created automatically)
│   ├── qc/                        # QC plots and filtered data
│   ├── de/                        # DE results and visualizations
│   ├── ora/                       # ORA results
│   └── gsea/                      # GSEA results
├── run_pipeline.R                 # Master script to run entire pipeline
└── README.md                      # This file
```

## Requirements

### R Version
- R >= 4.0.0

### Required R Packages

```r
# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "clusterProfiler",
  "ReactomePA",
  "org.Mm.eg.db",  # For mouse (change for other organisms)
  "apeglm",
  "ashr",
  "EnhancedVolcano",
  "pheatmap"
))

# Install from CRAN
install.packages(c(
  "tidyverse",
  "data.table",
  "ggrepel",
  "RColorBrewer",
  "gridExtra",
  "msigdbr"
))
```

### Organism-Specific Packages

For organisms other than mouse, install the appropriate annotation package:

```r
# Human
BiocManager::install("org.Hs.eg.db")

# Fly
BiocManager::install("org.Dm.eg.db")

# Rat
BiocManager::install("org.Rn.eg.db")
```

## Input Data Format

### 1. Gene Count Matrix (`data/gene_counts.tsv`)

Tab-separated file with genes as rows and samples as columns:

```
gene_id         gene_name    Sample1  Sample2  Sample3  Sample4
ENSMUSG00000000001  Gnai3       500      450      520      480
ENSMUSG00000000003  Pbsn        120      130      110      125
...
```

Requirements:
- First column: Gene IDs (Ensembl IDs recommended)
- Second column: Gene symbols/names
- Remaining columns: Raw read counts per sample
- Column names must match sample IDs in metadata file

### 2. Sample Metadata (`data/sample_metadata.tsv`)

Tab-separated file with sample information:

```
sample      condition
Sample1     Control
Sample2     Control
Sample3     Treatment1
Sample4     Treatment1
```

Requirements:
- Must include `sample` and `condition` columns
- Sample IDs must match count matrix column names
- Condition names used for differential expression contrasts

## Configuration

Edit `scripts/00_config.R` to customize analysis parameters:

```r
# Filtering thresholds
MIN_COUNT <- 10               # Minimum count per gene
MIN_SAMPLES <- 3              # Minimum samples with MIN_COUNT

# DE thresholds
PADJ_THRESHOLD <- 0.05        # Adjusted p-value cutoff
LFC_THRESHOLD <- 1            # Log2 fold change threshold (2-fold)

# Organism settings
ORGDB <- "org.Mm.eg.db"       # Organism database
ORGANISM <- "mouse"           # For ReactomePA
KEGG_ORGANISM <- "mmu"        # KEGG organism code
```

### Organism Codes

| Organism | ORGDB | ORGANISM | KEGG_ORGANISM | MSIGDB_SPECIES |
|----------|-------|----------|---------------|----------------|
| Mouse | org.Mm.eg.db | mouse | mmu | Mus musculus |
| Human | org.Hs.eg.db | human | hsa | Homo sapiens |
| Fly | org.Dm.eg.db | fly | dme | Drosophila melanogaster |
| Rat | org.Rn.eg.db | rat | rno | Rattus norvegicus |

## Usage

### Option 1: Run Complete Pipeline

Execute all analysis steps in sequence:

```bash
Rscript run_pipeline.R
```

### Option 2: Run Individual Steps

Run scripts separately for more control:

```r
# Step 1: QC and filtering
source("scripts/01_qc_filtering.R")

# Step 2: Differential expression
source("scripts/02_differential_expression.R")

# Step 3: Functional enrichment
source("scripts/03_functional_enrichment.R")
```

## Outputs

### Quality Control (`results/qc/`)

- `library_sizes.pdf` - Total read counts per sample
- `count_distribution.pdf` - Expression level distributions
- `dispersion.pdf` - Gene-wise dispersion estimates
- `PCA_prefiltering.pdf` - PCA before filtering
- `PCA_postfiltering.pdf` - PCA after filtering
- `sample_correlation.pdf` - Pearson correlation heatmap
- `filtered_data.RData` - Filtered DESeq2 object

### Differential Expression (`results/de/`)

- `*_all_genes.csv` - Results for all genes
- `*_sig_genes.csv` - Significantly differentially expressed genes
- `*_volcano.pdf` - Volcano plots
- `*_MA_plot.pdf` - MA plots
- `*_heatmap_top50.pdf` - Heatmaps of top DEGs
- `de_summary.csv` - Summary statistics table
- `de_results.RData` - Complete DE results

### ORA Analysis (`results/ora/`)

- `*_up_ORA_GO.csv` - GO enrichment for upregulated genes
- `*_down_ORA_GO.csv` - GO enrichment for downregulated genes
- `*_ORA_Reactome.csv` - Reactome pathway enrichment
- `*_ORA_KEGG.csv` - KEGG pathway enrichment
- `*_dotplot.pdf` - Dotplot visualizations

### GSEA Analysis (`results/gsea/`)

- `*_GSEA_GO.csv` - GO GSEA results
- `*_GSEA_Reactome.csv` - Reactome GSEA results
- `*_GSEA_KEGG.csv` - KEGG GSEA results
- `*_dotplot.pdf` - Dotplot visualizations
- `enrichment_results.RData` - Complete enrichment results

## Customization

### Defining Contrasts

Edit the contrasts in `scripts/02_differential_expression.R`:

```r
contrasts <- list(
  Treatment1_vs_Control = c("condition", "Treatment1", "Control"),
  Treatment2_vs_Control = c("condition", "Treatment2", "Control"),
  Treatment2_vs_Treatment1 = c("condition", "Treatment2", "Treatment1")
)
```

### Adjusting Color Scheme

Modify colors in `scripts/00_config.R`:

```r
ANN_COLORS <- list(
  condition = c(
    Control = "#D53E4F",      # Red
    Treatment1 = "#66C2A5",   # Teal
    Treatment2 = "#3288BD"    # Blue
  )
)
```

## Key Analysis Decisions

### Filtering Strategy

The pipeline uses `edgeR::filterByExpr()` for adaptive filtering, which:
- Keeps genes with sufficient reads in enough samples
- Adapts to experimental design and library sizes
- Generally removes ~50-70% of lowly-expressed genes

### Log Fold Change Shrinkage

The pipeline uses `ashr` for LFC shrinkage, which:
- Works with any contrast specification
- Reduces noise in lowly-expressed genes
- Improves ranking and visualization
- Alternative: `apeglm` (requires coefficient specification)

### GSEA Ranking Metric

Genes are ranked by: `sign(log2FC) * -log10(padj)`

This metric:
- Combines effect size (log2FC) and significance (padj)
- Preserves directionality (sign)
- Weights by statistical confidence
- Alternative: log2FC only (less recommended)

## Troubleshooting

### Issue: Too few/many significant genes

**Solution**: Adjust thresholds in `00_config.R`
```r
PADJ_THRESHOLD <- 0.1      # More lenient
LFC_THRESHOLD <- 0.585     # 1.5-fold change instead of 2-fold
```

### Issue: Enrichment analysis fails

**Possible causes**:
1. Too few significant genes → Lower DE thresholds
2. Gene ID conversion issues → Check organism database
3. Network issues → Try again or use cached databases

### Issue: Memory errors with large datasets

**Solution**: Process contrasts individually or increase R memory limit
```r
# Increase memory (Linux/Mac)
options(java.parameters = "-Xmx8g")
```

## Citation

If you use this pipeline, please cite:

- **DESeq2**: Love, MI et al. (2014) Genome Biology
- **clusterProfiler**: Yu G et al. (2012) OMICS
- **ReactomePA**: Yu G and He QY (2016) Molecular BioSystems

## Author

Tero Sievänen  
University of Eastern Finland  
Bioinformatics Center

## License

This pipeline is provided as-is for academic and research use.

## Version History

- v1.0 (2025-01): Initial modular release
