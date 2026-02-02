#!/usr/bin/env Rscript
################################################################################
# PCA Visualization Script
#
# Generates PCA plots from PLINK2 output:
# - Scree plot (variance explained per PC)
# - PC1 vs PC2 scatter plot
#
# Usage: Rscript plot_pca.R <eigenvec_file> <eigenval_file> <output_prefix>
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_pca.R <eigenvec_file> <eigenval_file> <output_prefix>", call. = FALSE)
}

eigenvec_file <- args[1]  
eigenval_file <- args[2]  
output_prefix <- args[3]  

cat("Loading PCA results...\n")

# Import PCA data
pca <- read_table(eigenvec_file, col_names = FALSE, show_col_types = FALSE)
eigenval <- scan(eigenval_file, quiet = TRUE)

# Clean and format data
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))
pca <- as_tibble(pca) %>%
  mutate(across(starts_with("PC"), as.numeric))

# Calculate percent variance explained (PVE)
pve <- data.frame(
  PC = 1:length(eigenval), 
  pve = eigenval / sum(eigenval) * 100
)

# Scree plot
cat("Generating scree plot...\n")
scree_plot_file <- paste0(output_prefix, "_scree.pdf")

p_scree <- ggplot(pve, aes(x = PC, y = pve)) +
  geom_point(size = 3, color = "steelblue") +  
  geom_line(group = 1, color = "steelblue", linewidth = 1) +  
  labs(
    title = "Scree Plot: Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

ggsave(scree_plot_file, p_scree, width = 10, height = 6, dpi = 300)

# PCA scatter plot (PC1 vs PC2)
cat("Generating PCA scatter plot...\n")
pca_plot_file <- paste0(output_prefix, "_scatter.pdf")

p_scatter <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, alpha = 0.6, color = "dodgerblue3") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +  
  coord_fixed(ratio = 1) +  
  labs(
    title = "Principal Component Analysis",
    x = sprintf("PC1 (%.2f%% variance)", pve$pve[1]),
    y = sprintf("PC2 (%.2f%% variance)", pve$pve[2])
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

ggsave(pca_plot_file, p_scatter, width = 8, height = 8, dpi = 300)

# Print cumulative variance
cumulative_pve <- cumsum(pve$pve)
cat("\nCumulative Variance Explained by First 10 PCs:\n")
print(head(cumulative_pve, 10))

cat("\nPCA plots saved:\n")
cat("  -", scree_plot_file, "\n")
cat("  -", pca_plot_file, "\n")
