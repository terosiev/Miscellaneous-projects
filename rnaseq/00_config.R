################################################################################
# Configuration File for RNA-seq Analysis Pipeline
#
# Define all paths, parameters, and settings here
# Source this file at the beginning of each analysis script
################################################################################

## PATHS ----

# Input data paths
DATA_DIR <- "data"
COUNTS_FILE <- file.path(DATA_DIR, "gene_counts.tsv")
METADATA_FILE <- file.path(DATA_DIR, "sample_metadata.tsv")

# Output directories
RESULTS_DIR <- "results"
QC_DIR <- file.path(RESULTS_DIR, "qc")
DE_DIR <- file.path(RESULTS_DIR, "de")
ORA_DIR <- file.path(RESULTS_DIR, "ora")
GSEA_DIR <- file.path(RESULTS_DIR, "gsea")

# Create output directories
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ORA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GSEA_DIR, recursive = TRUE, showWarnings = FALSE)

## ANALYSIS PARAMETERS ----

# Filtering thresholds
MIN_COUNT <- 10           # Minimum count per gene
MIN_SAMPLES <- 3          # Minimum number of samples with MIN_COUNT

# Differential expression thresholds
PADJ_THRESHOLD <- 0.05    # Adjusted p-value cutoff (FDR)
LFC_THRESHOLD <- 1        # Log2 fold change threshold (2-fold change)

# Enrichment analysis parameters
MIN_GENESET_SIZE <- 10    # Minimum genes in pathway
MAX_GENESET_SIZE <- 500   # Maximum genes in pathway
ENRICHMENT_PADJ <- 0.05   # Enrichment p-value cutoff

## ORGANISM SETTINGS ----

# Change these for different organisms:
# Mouse: org.Mm.eg.db, "mouse", "mmu"
# Human: org.Hs.eg.db, "human", "hsa"
# Fly: org.Dm.eg.db, "fly", "dme"

ORGDB <- "org.Mm.eg.db"           # Organism annotation database
ORGANISM <- "mouse"                # Organism name (for ReactomePA)
KEGG_ORGANISM <- "mmu"            # KEGG organism code
MSIGDB_SPECIES <- "Mus musculus"  # MSigDB species name

## COLOR PALETTE ----

# Define colors for experimental conditions
# Adjust these to match your experimental groups
ANN_COLORS <- list(
  condition = c(
    Control = "#D53E4F", 
    Treatment1 = "#66C2A5",
    Treatment2 = "#3288BD"
  )
)

## LOAD REQUIRED LIBRARIES ----

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(edgeR)
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(enrichplot)
})

# Load organism-specific database
if (ORGDB == "org.Mm.eg.db") {
  library(org.Mm.eg.db)
} else if (ORGDB == "org.Hs.eg.db") {
  library(org.Hs.eg.db)
} else if (ORGDB == "org.Dm.eg.db") {
  library(org.Dm.eg.db)
}

## SESSION INFO ----

message("Configuration loaded successfully")
message("Organism: ", ORGANISM)
message("Output directory: ", RESULTS_DIR)
