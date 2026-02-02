#!/usr/bin/env Rscript
################################################################################
# Genomic Inflation Factor (Lambda) Calculator
#
# Calculates lambda (genomic inflation factor) from GWAS results
# Lambda quantifies the extent of inflation in test statistics
#
# Usage: Rscript calculate_lambda.R <results_file> <output_file>
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript calculate_lambda.R <results_file> <output_file>")
}

results_file <- args[1]
output_file <- args[2]

# Check if file exists
if (!file.exists(results_file)) {
  stop("Results file not found: ", results_file)
}

cat("Reading GWAS results...\n")
results <- fread(results_file)

# Remove NAs
results <- results[!is.na(P)]

# Calculate lambda (all SNPs)
lambda_all <- median(qchisq(1 - results$P, 1), na.rm = TRUE) / 0.455

# Calculate lambda (non-significant SNPs only)
lambda_nonsig <- median(qchisq(1 - results$P[results$P > 5e-8], 1), na.rm = TRUE) / 0.455

# Count variants
n_total <- nrow(results)
n_sig <- sum(results$P < 5e-8, na.rm = TRUE)

# Print results
cat("\nGenomic Inflation Factor (Lambda):\n")
cat("  - All SNPs:", round(lambda_all, 4), "\n")
cat("  - Non-significant SNPs:", round(lambda_nonsig, 4), "\n")
cat("\nTotal variants:", n_total, "\n")
cat("Significant variants:", n_sig, "\n")

# Save lambda to file
write(round(lambda_all, 4), output_file)

cat("\nLambda saved to:", output_file, "\n")
