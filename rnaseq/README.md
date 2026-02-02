# Longitudinal RNA-seq Analysis Workflow

**Author:** terosiev
**Description:** Complete workflow for longitudinal RNA-seq differential expression analysis with quality control, pathway enrichment, and gene set enrichment analysis.

## Overview

This workflow implements a comprehensive analysis pipeline for longitudinal RNA-seq data using limma-voom with duplicate correlation to account for repeated measures on the same subjects. The pipeline includes:

1. **Quality Control** (`qc.R`) - Data preprocessing, filtering, and normalization
2. **Differential Expression** (`de.R`) - Longitudinal DE analysis with limma-voom
3. **Over-Representation Analysis** (`ora.R`) - Gene Ontology enrichment analysis
4. **Gene Set Enrichment Analysis** (`gsea.R`) - GSEA using MSigDB gene sets

## Key Features

- **Longitudinal Design**: Handles repeated measures using duplicate correlation
- **Batch Correction**: Includes batch effects in the model design
- **Comprehensive QC**: PCA plots, hierarchical clustering, and sample correlations
- **Multiple Contrasts**: Supports complex experimental designs with interactions
- **Rich Visualizations**: Volcano plots, heatmaps, dotplots, cnetplots, and enrichment plots

## Requirements

### R Packages

```r
# Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "edgeR", 
  "limma",
  "tximport",
  "clusterProfiler",
  "enrichplot",
  "AnnotationHub",
  "ensembldb",
  "org.Hs.eg.db",
  "DOSE",
  "pathview"
))

# CRAN packages
install.packages(c(
  "tidyverse",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "EnhancedVolcano",
  "cowplot",
  "msigdbr",
  "ggplotify"
))
```

## Input Data Structure

The workflow expects the following input files:

```
project/
├── data/
│   ├── counts_matrix.csv      # Gene count matrix (genes x samples)
│   ├── metadata.csv            # Sample metadata with experimental design
│   └── tx2gene.tsv             # Transcript-to-gene mapping (GrCh38)
└── results/
    ├── qc/
    ├── differential_expression/
    ├── ora/
    └── gsea/
```

### Metadata Format

Your `metadata.csv` should contain:
- `sample_id` - Unique sample identifier
- `subject_id` - Subject/individual identifier (for repeated measures)
- `timepoint` - Time point (e.g., time1, time2, time3)
- `condition` - Experimental condition
- `batch` - Sequencing batch

## Usage

### 1. Quality Control

```r
source("qc.R")
```

**Outputs:**
- PCA plots colored by batch, time, and condition
- Sample correlation heatmap
- Expression density plots (before/after filtering)
- MDS plots

### 2. Differential Expression Analysis

```r
source("de.R")
```

**Outputs:**
- DEG lists (significant genes, adj.P < 0.05)
- Complete results for all genes
- Volcano plots for each contrast
- Heatmaps of top DEGs (up to 100 genes)

### 3. Over-Representation Analysis

```r
source("ora.R")
```

**Outputs:**
- GO enrichment results (all DEGs, upregulated, downregulated)
- Dotplots showing enriched terms
- Cnetplots connecting genes to pathways

### 4. Gene Set Enrichment Analysis

```r
source("gsea.R")
```

**Outputs:**
- GSEA results using MSigDB gene sets
- Dotplots of enriched gene sets
- GSEA enrichment plots (top 10)
- Ridgeplots showing enrichment distributions
- Cnetplots for top up/downregulated pathways

## Customization

### Modifying Contrasts

Edit the contrast matrix in `de_anonymized.R` to match your experimental design:

```r
cm <- makeContrasts(
  # Within-condition time comparisons
  Cond1_T2T1 = groupCONDITION1_time2 - groupCONDITION1_time1,
  
  # Between-condition comparisons
  Cond1_Cond2_T1 = groupCONDITION1_time1 - groupCONDITION2_time1,
  
  # Interaction terms
  Interaction_T2T1 = (groupCONDITION1_time2 - groupCONDITION1_time1) - 
                     (groupCONDITION2_time2 - groupCONDITION2_time1),
  
  levels=design
)
```

### Filtering Thresholds

Adjust gene filtering parameters in `qc_anonymized.R`:

```r
keep.exprs <- edgeR::filterByExpr(data, 
                                   min.count = 30, 
                                   min.total.count = 50)
```

### Significance Cutoffs

Modify thresholds in respective scripts:
- **DE analysis**: `p = 0.05`, `logFC > 1`
- **ORA**: `qvalueCutoff = 0.01`, `minGSSize = 20`
- **GSEA**: `pvalueCutoff = 0.05`, `minGSSize = 10`

## Statistical Methods

### Duplicate Correlation

For longitudinal data with repeated measures, the workflow uses `duplicateCorrelation()` from limma to estimate within-subject correlation:

```r
corfit <- duplicateCorrelation(v, design, block=data$samples$subject)
```

This correlation is then incorporated into the linear model:

```r
fit <- lmFit(v, design, block=data$samples$subject, 
             correlation=corfit$consensus)
```

### Normalization

- TMM normalization (trimmed mean of M-values) for count data
- voom transformation for variance stabilization
- rlog transformation (DESeq2) for visualization

## Tips

1. **Sample Size**: Duplicate correlation works best with 3+ replicates per subject
2. **Batch Effects**: Always include batch in your design if samples were processed in multiple batches
3. **Filtering**: The `filterByExpr()` function automatically adjusts based on library sizes and design
4. **Memory**: GSEA with large gene sets may require substantial RAM (8GB+ recommended)

## Citation

If you use this workflow, please cite the key packages:

- **limma**: Ritchie ME, et al. (2015). Nucleic Acids Research.
- **edgeR**: Robinson MD, et al. (2010). Bioinformatics.
- **clusterProfiler**: Yu G, et al. (2012). OMICS.

## License

This workflow is provided under the MIT License. Feel free to modify and adapt for your own research.
