#!/usr/bin/env Rscript
################################################################################
# GWAS Visualization Script
#
# Generates Manhattan and QQ plots from GWAS results
#
# Usage: Rscript plot_gwas.R <results_file> <output_dir> <plot_prefix>
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(qqman)
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_gwas.R <results_file> <output_dir> <plot_prefix>")
}

results_file <- args[1]
output_dir <- args[2]
plot_prefix <- args[3]

# Check if file exists
if (!file.exists(results_file)) {
  stop("Results file not found: ", results_file)
}

cat("Reading GWAS results:", results_file, "\n")
results <- fread(results_file)

# Filter to autosomes only
results <- results[`#CHROM` %in% 1:22]

cat("Total variants:", nrow(results), "\n")

# Count significant variants
n_sig <- sum(results$P < 5e-8, na.rm = TRUE)
n_sug <- sum(results$P >= 5e-8 & results$P < 1e-5, na.rm = TRUE)

cat("Genome-wide significant (P < 5e-8):", n_sig, "\n")
cat("Suggestive (P < 1e-5):", n_sug, "\n")

# Manhattan plot
cat("\nGenerating Manhattan plot...\n")
manhattan_file <- file.path(output_dir, paste0("manhattan_", plot_prefix, ".pdf"))

pdf(manhattan_file, width = 14, height = 7)
manhattan(results,
          chr = "#CHROM",
          bp = "POS",
          p = "P",
          snp = "ID",
          main = paste("Manhattan Plot:", plot_prefix),
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5),
          col = c("blue4", "orange3"),
          chrlabs = as.character(1:22),
          annotatePval = 5e-8,
          annotateTop = TRUE)
dev.off()

# QQ plot
cat("Generating QQ plot...\n")
qq_file <- file.path(output_dir, paste0("qq_", plot_prefix, ".pdf"))

pdf(qq_file, width = 7, height = 7)
qq(results$P,
   main = paste("QQ Plot:", plot_prefix))
dev.off()

cat("\nPlots saved:\n")
cat("  - Manhattan:", manhattan_file, "\n")
cat("  - QQ plot:", qq_file, "\n")
