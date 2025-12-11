################################################################################
# RNA-seq Analysis Pipeline
# 
# Complete workflow for RNA-seq differential expression analysis including:
# 1. Quality control and data filtering
# 2. Differential expression analysis with DESeq2
# 3. Functional enrichment analysis (ORA and GSEA)
#
# Assumes preprocessing with nf-core/rnaseq or similar pipeline
################################################################################


################################################################################
# SCRIPT 1: Quality Control and Data Filtering
################################################################################

# Load required libraries
library(data.table)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(edgeR)

## SETUP ----

# Define output directory
output_path_qc <- "results/qc/"
dir.create(output_path_qc, recursive = TRUE, showWarnings = FALSE)

# Define color palette for experimental conditions
# Adjust these to match your experimental groups
ann_colors <- list(
  condition = c(
    Control = "#D53E4F", 
    Treatment1 = "#66C2A5",
    Treatment2 = "#3288BD"
  )
)

## DATA IMPORT ----

# Import raw count matrix from salmon/STAR/etc
# Rows = genes, columns = samples
counts.df <- data.frame(fread("data/gene_counts.tsv", sep = "\t"))

# Import sample metadata
# Must include: sample, condition columns
meta <- data.frame(fread("data/sample_metadata.tsv", sep = "\t"))

# Store gene annotation (gene IDs and symbols)
gene_annotation <- counts.df %>% 
  dplyr::select(gene_id)
gene_annotation$symbol <- counts.df$gene_name  # Or use gene_id if symbols not available

# Prepare count matrix
# Set gene IDs as rownames and remove annotation columns
rownames(counts.df) <- counts.df$gene_id
counts.df <- counts.df %>% 
  dplyr::select(-gene_id, -gene_name)  # Remove any non-count columns

# Round counts to integers (required for DESeq2)
counts.df <- round(counts.df)

# Prepare metadata
# Set sample IDs as rownames
rownames(meta) <- meta$sample

## DATA INSPECTION ----

# Summary statistics of raw counts
summary(counts.df)

# Calculate library sizes (total reads per sample)
library_sizes <- data.frame(
  Sample = colnames(counts.df),
  TotalCounts = colSums(counts.df),
  Condition = meta$condition[match(colnames(counts.df), rownames(meta))]
)

# Calculate mean and median library sizes
mean_counts <- mean(library_sizes$TotalCounts)
median_counts <- median(library_sizes$TotalCounts)

# Visualize library sizes
# Helps identify samples with unusually low/high sequencing depth
p <- ggplot(library_sizes, aes(x = Sample, y = TotalCounts, fill = Condition)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = mean_counts, linetype = "Mean"), 
             color = "black", linewidth = 1) +
  geom_hline(aes(yintercept = median_counts, linetype = "Median"), 
             color = "gray40", linewidth = 1) +
  scale_fill_manual(values = ann_colors$condition) +
  scale_linetype_manual(name = "Reference", 
                       values = c("Mean" = "solid", "Median" = "dashed")) +
  ggtitle("Library Sizes") +
  xlab("Sample") +
  ylab("Total Counts") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  annotate("text", x = 1, y = mean_counts * 1.05, 
           label = sprintf("Mean: %.0f", mean_counts), hjust = 0) +
  annotate("text", x = 1, y = median_counts * 0.95, 
           label = sprintf("Median: %.0f", median_counts), hjust = 0)

ggsave(file.path(output_path_qc, "library_sizes.pdf"), plot = p, 
       width = 8, height = 6)

# Visualize count distribution across samples
# Log2 transformation for better visualization
count_dist <- counts.df %>%
  mutate(Gene = rownames(counts.df)) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Counts") %>%
  mutate(LogCounts = log2(Counts + 1)) %>%  # +1 to avoid log(0)
  left_join(dplyr::select(meta, sample, condition), by = c("Sample" = "sample"))

# Calculate mean and median log counts
mean_logcounts <- mean(count_dist$LogCounts)
median_logcounts <- median(count_dist$LogCounts)

# Plot count distribution
# Boxplots show median, quartiles, and outliers per sample
p <- ggplot(count_dist, aes(x = Sample, y = LogCounts, fill = condition)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = mean_logcounts, linetype = "Mean"), 
             color = "black", linewidth = 1) +
  geom_hline(aes(yintercept = median_logcounts, linetype = "Median"), 
             color = "gray40", linewidth = 1) +
  scale_fill_manual(values = ann_colors$condition) +
  scale_linetype_manual(name = "Reference", 
                       values = c("Mean" = "solid", "Median" = "dashed")) +
  ggtitle("Log2 Count Distribution") +
  xlab("Sample") +
  ylab("Log2(Counts + 1)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave(file.path(output_path_qc, "count_distribution.pdf"), plot = p, 
       width = 8, height = 6)

## CREATE DESEQ2 OBJECT ----

# Prepare column data (sample metadata)
colData <- meta %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "sample") %>%
  dplyr::select(condition)

# Create DESeq2 dataset
# Design formula: ~ condition (tests for differences between conditions)
dds <- DESeqDataSetFromMatrix(
  countData = counts.df, 
  colData = colData, 
  design = ~ condition
)

# Estimate size factors (normalization for sequencing depth)
dds <- estimateSizeFactors(dds)

# Estimate dispersion (gene-wise variability)
# Dispersion = variance / mean, measures biological variability
dds <- estimateDispersions(dds)

# Plot dispersion estimates
# Black dots = gene-wise estimates
# Red line = fitted trend
# Blue circles = final estimates after shrinkage
pdf(file.path(output_path_qc, "dispersion.pdf"))
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

## PCA ANALYSIS (PRE-FILTERING) ----

# Function to perform and plot PCA
# PCA reveals sample clustering and potential batch effects
run_PCA <- function(dds, design_matrix, pdf_filename = NULL, path = NULL) {
  
  # Variance stabilizing transformation (VST)
  # Transforms counts to log2-like scale with stabilized variance
  # blind = TRUE: ignore design matrix (for exploratory analysis)
  rld <- vst(dds, blind = TRUE)
  
  # Perform PCA on transposed data (samples as rows)
  pca <- prcomp(t(assay(rld)))
  
  # Calculate variance explained by each PC
  var_explained <- summary(pca)$importance[2, ] * 100
  
  # Combine metadata with PC coordinates
  df <- data.frame(cbind(design_matrix, pca$x))
  
  # Define grouping variables for coloring
  intgroup <- c("condition")
  
  # Create PCA plots for each grouping variable
  plots <- lapply(intgroup, function(group) {
    plot <- ggplot(df, aes(x = PC1, y = PC2, color = .data[[group]])) + 
      geom_point(size = 7, alpha = 0.8) +
      geom_point(size = 7, shape = 1, color = "black") +  # Black outline
      geom_text(aes(label = rownames(design_matrix)), 
               nudge_x = 0.8, nudge_y = 0.8, 
               check_overlap = TRUE, color = "black") +
      theme_classic() +
      ggtitle(paste("PCA Plot (Group:", group, ")")) +
      xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
      ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      scale_color_manual(values = ann_colors$condition)
    
    # Save plot if filename provided
    if (!is.null(pdf_filename)) {
      filename <- if (!is.null(path)) {
        file.path(path, paste0(pdf_filename, "_", group, ".pdf"))
      } else {
        paste0(pdf_filename, "_", group, ".pdf")
      }
      tryCatch({
        ggsave(filename, plot, dpi = 300)
      }, error = function(e) {
        message(paste("Error saving plot:", e$message))
      })
    }
    
    return(plot)
  })
  
  return(plots)
}

# Run PCA before filtering
pdf_filename <- "PCA_prior_filtering"
run_PCA(dds, meta, pdf_filename, output_path_qc)

## SAMPLE CORRELATION ANALYSIS (PRE-FILTERING) ----

# Compute VST-transformed counts
rld <- vst(dds, blind = TRUE)
rld_mat <- assay(rld)

# Calculate pairwise sample correlations
# High correlations expected between biological replicates
rld_cor <- cor(rld_mat)

# Summary of correlation values
summary(as.vector(rld_cor[upper.tri(rld_cor)]))
print(round(rld_cor, 3))

# Prepare annotation for heatmap
anno <- data.frame(condition = meta$condition, row.names = rownames(meta))

# Define color palette for heatmap (red-white-blue)
heat_colors <- colorRampPalette(brewer.pal(9, "RdBu"))(100)

# Heatmap of sample correlations
# Samples with similar expression profiles cluster together
pheatmap(rld_cor, 
         annotation = anno, 
         annotation_colors = ann_colors, 
         cluster_rows = FALSE,  # Keep sample order
         cluster_cols = FALSE, 
         angle_col = 45, 
         fontsize = 7,
         filename = file.path(output_path_qc, "qc_hmap_prior.pdf"))

# Sample distance heatmap (Euclidean distance)
# Alternative to correlation; shows overall similarity
sample_dist <- dist(t(rld_mat))
pheatmap(as.matrix(sample_dist),
         cluster_cols = TRUE,  # Hierarchical clustering
         annotation = anno, 
         annotation_colors = ann_colors, 
         filename = file.path(output_path_qc, "qc_dist_prior.pdf"))

# Correlation using top 1000 most variable genes
# Focuses on genes driving most variation
rv <- rowVars(rld_mat)
top_genes <- order(rv, decreasing = TRUE)[1:1000]
rld_cor_top <- cor(rld_mat[top_genes, ])
pheatmap(rld_cor_top, 
         annotation = anno, 
         cluster_cols = FALSE,
         annotation_colors = ann_colors, 
         filename = file.path(output_path_qc, "qc_hmap_top1000_prior.pdf"))

## DATA FILTERING ----

# Convert to matrix and validate
counts <- as.matrix(counts.df)

# Check for duplicate or missing gene names
if (anyDuplicated(rownames(counts)) > 0 || any(is.na(rownames(counts)))) {
  stop("Row names of counts matrix are not unique or contain NAs")
}

# Verify metadata alignment with count matrix
if (!all(rownames(meta) == colnames(counts))) {
  stop("Meta row names do not match counts column names")
}

# Report pre-filtering gene count
cat("Genes before filtering:", nrow(dds), "\n")

# Filter low-expression genes
# Removes genes with insufficient counts across samples
# filterByExpr: adaptive filtering based on library sizes and group sizes
# min.count: minimum count threshold per sample
# min.total.count: minimum total count across all samples
# large.n: minimum number of samples with counts >= min.count
keep <- filterByExpr(
  counts(dds), 
  group = meta$condition, 
  min.count = 50,        # At least 50 reads in individual samples
  min.total.count = 100, # At least 100 reads total
  large.n = 3            # In at least 3 samples
)

# Apply filter
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# Add gene symbols to DESeq object
mcols(dds)$symbol <- gene_annotation$symbol[match(rownames(dds), 
                                                   gene_annotation$gene_id)]

# Calculate percentage of genes retained
pct_retained <- nrow(dds) / nrow(counts) * 100
cat("Percentage of genes retained:", round(pct_retained, 2), "%\n")

## VISUALIZE FILTERING EFFECT ----

# Compute log-CPM (counts per million) for visualization
# CPM normalizes for library size differences
lcpm_raw <- cpm(counts, log = TRUE, prior.count = 2)
lcpm_filt <- cpm(counts(dds), log = TRUE, prior.count = 2)

# Calculate filtering cutoff line
L <- mean(colSums(counts)) * 1e-6  # Average library size in millions
M <- median(colSums(counts)) * 1e-6  # Median library size
lcpm.cutoff <- log2(10 / M + 2 / L)

# Prepare data for density plots
df_raw <- data.frame(
  logCPM = as.vector(lcpm_raw),
  Sample = rep(colnames(counts), each = nrow(counts)),
  Condition = rep(meta$condition, each = nrow(counts))
)
df_filt <- data.frame(
  logCPM = as.vector(lcpm_filt),
  Sample = rep(colnames(counts), each = nrow(dds)),
  Condition = rep(meta$condition, each = nrow(dds))
)

# Density plot: raw data
# Shows distribution of gene expression levels
p1 <- ggplot(df_raw, aes(x = logCPM, color = Condition)) +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = lcpm.cutoff, linetype = "dashed", color = "black") +
  scale_color_manual(values = ann_colors$condition) +
  ggtitle("A. Raw Data") +
  xlab("Log-CPM") +
  theme_classic() +
  theme(legend.position = "top")

# Density plot: filtered data
# Removal of low-expression genes shifts distribution right
p2 <- ggplot(df_filt, aes(x = logCPM, color = Condition)) +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = lcpm.cutoff, linetype = "dashed", color = "black") +
  scale_color_manual(values = ann_colors$condition) +
  ggtitle("B. Filtered Data") +
  xlab("Log-CPM") +
  theme_classic() +
  theme(legend.position = "top")

# Ensure consistent y-axis limits
max_y <- max(c(
  sapply(1:ncol(lcpm_raw), function(i) max(density(lcpm_raw[, i])$y)),
  sapply(1:ncol(lcpm_filt), function(i) max(density(lcpm_filt[, i])$y))
))
p1 <- p1 + coord_cartesian(ylim = c(0, max_y * 1.1))
p2 <- p2 + coord_cartesian(ylim = c(0, max_y * 1.1))

# Combine and save density plots
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave(file.path(output_path_qc, "filtering.pdf"), plot = combined_plot, 
       width = 10, height = 5)

## PCA AND CORRELATIONS (POST-FILTERING) ----

# PCA after filtering
pdf_filename <- "PCA_post_filtering"
run_PCA(dds, meta, pdf_filename, output_path_qc)

# Sample correlations after filtering
rld <- vst(dds, blind = FALSE)  # Use design for better variance estimates
rld_cor <- cor(assay(rld))
pheatmap(rld_cor, 
         annotation = anno, 
         annotation_colors = ann_colors, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         angle_col = 45, 
         fontsize = 7, 
         filename = file.path(output_path_qc, "qc_hmap_post.pdf"))

# Sample distances after filtering
sample_dist <- dist(t(assay(vst(dds, blind = FALSE))))
pheatmap(as.matrix(sample_dist), 
         annotation = anno, 
         cluster_cols = TRUE,
         annotation_colors = ann_colors, 
         filename = file.path(output_path_qc, "qc_dist_post.pdf"))

# Top variable genes after filtering
rld_mat_post <- assay(rld)
rv <- rowVars(rld_mat_post)
top_genes <- order(rv, decreasing = TRUE)[1:1000]
rld_cor_top <- cor(rld_mat_post[top_genes, ])
pheatmap(rld_cor_top, 
         annotation = anno, 
         cluster_cols = FALSE,
         annotation_colors = ann_colors,
         filename = file.path(output_path_qc, "qc_hmap_top1000_post.pdf"))

## VISUALIZE NORMALIZATION EFFECT ----

# Compute VST without normalization (size factors = 1)
sizeFactors(dds) <- rep(1, ncol(dds))
vst_unnorm <- vst(dds, blind = TRUE)

# Compute VST with proper normalization
dds <- estimateSizeFactors(dds)
vst_norm <- vst(dds, blind = TRUE)

# Prepare data for boxplots
df_unnorm <- data.frame(
  logCounts = as.vector(assay(vst_unnorm)),
  Sample = rep(colnames(dds), each = nrow(dds)),
  Condition = rep(meta$condition, each = nrow(dds))
)
df_norm <- data.frame(
  logCounts = as.vector(assay(vst_norm)),
  Sample = rep(colnames(dds), each = nrow(dds)),
  Condition = rep(meta$condition, each = nrow(dds))
)

# Boxplot: unnormalized data
# Shows raw differences in library sizes
p3 <- ggplot(df_unnorm, aes(x = Sample, y = logCounts, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = ann_colors$condition) +
  ggtitle("A. Unnormalized Data") +
  ylab("VST Counts") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot: normalized data
# Medians should align after normalization
p4 <- ggplot(df_norm, aes(x = Sample, y = logCounts, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = ann_colors$condition) +
  ggtitle("B. Normalized Data") +
  ylab("VST Counts") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine and save boxplots
combined_boxplot <- grid.arrange(p3, p4, ncol = 2)
ggsave(file.path(output_path_qc, "normalization.pdf"), plot = combined_boxplot, 
       width = 10, height = 5)

# Display normalization factors (size factors)
cat("\nSize factors:\n")
print(sizeFactors(dds))


################################################################################
# SCRIPT 2: Differential Expression Analysis
################################################################################

# Load required libraries
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(ashr)

## SETUP ----

# Define output directory
output_path_deg <- "results/de/"
dir.create(output_path_deg, recursive = TRUE, showWarnings = FALSE)

# Use same color palette as QC
ann_colors <- list(
  condition = c(
    Control = "#D53E4F", 
    Treatment1 = "#66C2A5",
    Treatment2 = "#3288BD"
  )
)

## RUN DESEQ2 ----

# Perform differential expression analysis
# This runs the full DESeq2 pipeline:
# 1. Estimate size factors (normalization)
# 2. Estimate dispersions (variance)
# 3. Fit negative binomial GLM
# 4. Perform Wald test for each coefficient
dds <- DESeq(dds)

# Check available contrasts
print("Available coefficients in resultsNames(dds):")
print(resultsNames(dds))

## DEFINE CONTRASTS ----

# Define pairwise comparisons to test
# Each contrast compares two conditions
contrasts <- list(
  Treatment1_vs_Control = list(
    contrast = c("condition", "Treatment1", "Control"),
    coef = "condition_Treatment1_vs_Control"  # Must match resultsNames(dds)
  ),
  Treatment2_vs_Control = list(
    contrast = c("condition", "Treatment2", "Control"),
    coef = "condition_Treatment2_vs_Control"
  ),
  Treatment2_vs_Treatment1 = list(
    contrast = c("condition", "Treatment2", "Treatment1"),
    coef = "condition_Treatment2_vs_Treatment1"
  )
)

## PERFORM DE ANALYSIS FOR EACH CONTRAST ----

# Initialize storage lists
res_shrunken_sig_list <- list()  # Significant genes only
res_shrunken_all_list <- list()  # All genes
res_shrunken_list <- list()      # DESeq results objects

# Initialize summary statistics table
de_summary <- data.frame(
  Contrast = character(),
  Total_Significant = integer(),
  Up_Regulated = integer(),
  Down_Regulated = integer(),
  stringsAsFactors = FALSE
)

# Loop through each contrast
for (contrast_name in names(contrasts)) {
  cat("\n=== Processing contrast:", contrast_name, "===\n")
  
  # Extract contrast and coefficient name
  contrast <- contrasts[[contrast_name]]$contrast
  coef <- contrasts[[contrast_name]]$coef
  
  # Get unshrunken results
  # alpha: FDR cutoff for filtering
  res <- results(dds, contrast = contrast, alpha = 0.05)
  
  # Shrink log2 fold changes
  # Shrinkage reduces noise in lowly-expressed genes
  # Two methods available:
  # - apeglm: adaptive shrinkage (recommended, requires coef)
  # - ashr: adaptive shrinkage (works with contrasts)
  res_shrunken <- NULL
  if (!is.null(coef) && coef %in% resultsNames(dds)) {
    # Use apeglm if coefficient exists
    res_shrunken <- lfcShrink(dds, coef = coef, type = "apeglm", res = res)
  } else {
    # Use ashr as fallback
    message(paste("Coefficient for", contrast_name, 
                 "not found. Using ashr for lfcShrink."))
    res_shrunken <- lfcShrink(dds, contrast = contrast, type = "ashr", res = res)
  }
  
  # Add gene symbols
  res_shrunken$symbol <- rownames(res_shrunken)
  
  # Convert to data frame
  res_shrunken_all <- as.data.frame(res_shrunken)
  res_shrunken_all$gene <- rownames(res_shrunken_all)
  
  # Filter for significant genes
  # Significance criteria:
  # - Adjusted p-value < 0.05 (FDR control)
  # - |log2 fold change| > 1.0 (at least 2-fold change)
  res_shrunken_sig <- res_shrunken_all %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1.0) %>%
    arrange(padj)
  
  # Calculate summary statistics
  total_sig <- nrow(res_shrunken_sig)
  up_reg <- sum(res_shrunken_sig$log2FoldChange > 0, na.rm = TRUE)
  down_reg <- sum(res_shrunken_sig$log2FoldChange < 0, na.rm = TRUE)
  
  cat("Total significant genes:", total_sig, "\n")
  cat("Up-regulated:", up_reg, "\n")
  cat("Down-regulated:", down_reg, "\n")
  
  # Append to summary table
  de_summary <- rbind(de_summary, data.frame(
    Contrast = contrast_name,
    Total_Significant = total_sig,
    Up_Regulated = up_reg,
    Down_Regulated = down_reg
  ))
  
  # Store results
  res_shrunken_sig_list[[contrast_name]] <- res_shrunken_sig
  res_shrunken_all_list[[contrast_name]] <- res_shrunken_all
  res_shrunken_list[[contrast_name]] <- res_shrunken
  
  # Save results to CSV
  write.csv(res_shrunken_all, 
           file.path(output_path_deg, paste0(contrast_name, "_all_genes.csv")), 
           row.names = FALSE)
  
  if (nrow(res_shrunken_sig) > 0) {
    write.csv(res_shrunken_sig, 
             file.path(output_path_deg, paste0(contrast_name, "_sig_genes.csv")), 
             row.names = FALSE)
  }
}

# Save summary statistics
write.csv(de_summary, file.path(output_path_deg, "de_summary_stats.csv"), 
         row.names = FALSE)

## VISUALIZE RESULTS ----

# 1. Log2 Fold Change Distribution Histograms
for (contrast_name in names(contrasts)) {
  p <- ggplot(res_shrunken_sig_list[[contrast_name]], 
             aes(x = log2FoldChange)) +
    geom_histogram(bins = 50, fill = "#3288BD", color = "black") +
    ggtitle(paste("Log2FC Distribution:", gsub("_", " ", contrast_name))) +
    xlab("Log2 Fold Change") +
    ylab("Count") +
    theme_classic()
  
  ggsave(file.path(output_path_deg, 
                  paste0(contrast_name, "_log2FC_histogram.pdf")), 
        plot = p, width = 6, height = 4, dpi = 300)
}

# 2. MA Plots (log2FC vs mean expression)
# Shows relationship between fold change and expression level
# Points in red: significant genes (padj < 0.05)
for (contrast_name in names(contrasts)) {
  res_ma <- res_shrunken_list[[contrast_name]]
  
  pdf(file.path(output_path_deg, paste0(contrast_name, "_MAplot.pdf")), 
      width = 8, height = 6)
  DESeq2::plotMA(res_ma, 
                main = paste("MA Plot:", gsub("_", " ", contrast_name)), 
                ylim = c(-5, 5),
                alpha = 0.05,
                colSig = "red",
                colNonSig = "black")
  dev.off()
}

# 3. Volcano Plots (log2FC vs -log10 p-value)
# Shows both magnitude and significance of changes
for (contrast_name in names(contrasts)) {
  contrast_data <- res_shrunken_all_list[[contrast_name]]
  
  # Calculate y-axis limit
  max_y <- max(-log10(contrast_data$padj[!is.na(contrast_data$padj)]), 
               na.rm = TRUE)
  ylim_max <- min(ceiling(max_y * 1.1), 20)  # Cap at 20
  
  # Calculate x-axis limit
  max_fc <- max(abs(contrast_data$log2FoldChange[
    !is.na(contrast_data$log2FoldChange)]), na.rm = TRUE)
  xlim_max <- ceiling(max_fc * 1.1)
  
  # Create volcano plot
  volcano_plot <- EnhancedVolcano(
    contrast_data,
    x = "log2FoldChange",
    y = "padj",
    lab = contrast_data$symbol,
    pCutoff = 0.05,           # Significance threshold
    FCcutoff = 1.0,           # Fold change threshold
    labSize = 3,
    colAlpha = 0.8,
    col = c("grey30", "#66C2A5", "#3288BD", "#D53E4F"),
    title = paste("DE Genes:", gsub("_", " ", contrast_name)),
    subtitle = NULL,
    legendPosition = "bottom",
    legendLabels = c("NS", "Log2 FC", "P-value", "Sig"),
    drawConnectors = FALSE
  ) +
    coord_cartesian(ylim = c(0, ylim_max), 
                   xlim = c(-xlim_max, xlim_max)) +
    theme_classic()
  
  ggsave(file.path(output_path_deg, paste0(contrast_name, "_volcano.pdf")), 
        plot = volcano_plot, width = 8, height = 6, dpi = 300)
}

# 4. Heatmaps of Significant Genes
# Shows expression patterns of top DE genes across samples

# Get normalized counts (VST-transformed)
vst_norm <- vst(dds, blind = FALSE)
vst_counts <- assay(vst_norm)
normalized_counts <- as.data.frame(vst_counts) %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(symbol = mcols(dds)$symbol)

# Create heatmaps for each contrast
for (contrast_name in names(contrasts)) {
  # Extract conditions for current contrast
  conditions <- contrasts[[contrast_name]]$contrast[2:3]
  
  # Filter normalized counts for relevant samples only
  sample_conditions <- meta$condition
  sample_names <- colnames(vst_counts)
  keep_samples <- sample_names[sample_conditions %in% conditions]
  norm_counts_filtered <- normalized_counts %>%
    dplyr::select(gene, symbol, all_of(keep_samples))
  
  # Filter metadata for relevant samples
  meta_subset <- meta[meta$condition %in% conditions, , drop = FALSE]
  anno <- data.frame(condition = meta_subset$condition, 
                    row.names = rownames(meta_subset))
  
  # Color palette for heatmap
  heat_colors <- rev(brewer.pal(7, "RdBu"))
  
  # Get significant genes for this contrast
  sig_genes <- res_shrunken_sig_list[[contrast_name]]$gene
  if (length(sig_genes) == 0) {
    message(paste("No significant genes for", contrast_name, "- skipping heatmap"))
    next
  }
  
  # Select top 100 most significant genes
  norm_sig <- norm_counts_filtered %>%
    dplyr::filter(gene %in% sig_genes) %>%
    dplyr::left_join(res_shrunken_sig_list[[contrast_name]] %>% 
                    dplyr::select(gene, padj), by = "gene") %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 100)
  
  if (nrow(norm_sig) == 0) {
    message(paste("No data for heatmap:", contrast_name))
    next
  }
  
  # Prepare matrix for heatmap
  heatmap_data <- norm_sig %>%
    column_to_rownames("gene") %>%
    dplyr::select(-padj, -symbol) %>%
    as.matrix()
  
  # Z-score transformation (row-wise scaling)
  # Centers and scales each gene across samples
  heatmap_zscore <- t(scale(t(heatmap_data)))
  
  # Subset annotation colors for current conditions
  ann_colors_subset <- list(
    condition = ann_colors$condition[names(ann_colors$condition) %in% conditions]
  )
  
  # Create heatmap
  pdf(file.path(output_path_deg, paste0(contrast_name, "_heatmap.pdf")), 
      width = 8, height = 8)
  pheatmap(heatmap_zscore,
          color = heat_colors,
          cluster_rows = TRUE,    # Cluster genes by expression pattern
          cluster_cols = FALSE,   # Keep sample order
          show_rownames = TRUE,   # Show gene names
          annotation_col = anno,
          annotation_colors = ann_colors_subset,
          fontsize = 5,
          fontsize_row = 5,
          angle_col = 45,
          main = paste("Top DE Genes:", gsub("_", " ", contrast_name)))
  dev.off()
}


################################################################################
# SCRIPT 3: Functional Enrichment Analysis (ORA and GSEA)
################################################################################

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse annotations (change for other organisms)
library(DOSE)
library(enrichplot)
library(tidyverse)
library(ReactomePA)
library(msigdbr)

## SETUP ----

# Define output directories
output_path_ora <- "results/ora/"
output_path_gsea <- "results/gsea/"
dir.create(output_path_ora, recursive = TRUE, showWarnings = FALSE)
dir.create(output_path_gsea, recursive = TRUE, showWarnings = FALSE)

## RETRIEVE GENE SETS ----

# Get all gene sets from MSigDB
# MSigDB: Molecular Signatures Database (curated gene sets)
all_gene_sets <- msigdbr(species = "Mus musculus")  # Change for other species

# Get KEGG pathways
# KEGG: Kyoto Encyclopedia of Genes and Genomes (metabolic/signaling pathways)
kegg_sets <- msigdbr(species = "Mus musculus", 
                    category = "C2", 
                    subcategory = "CP:KEGG_MEDICUS")

# Format KEGG gene sets for clusterProfiler
kegg_gmt <- kegg_sets %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(entrez_gene)) %>%
  dplyr::mutate(genes = lapply(genes, as.character)) %>%
  dplyr::select(gs_name, genes)

## OVERREPRESENTATION ANALYSIS (ORA) ----

# ORA tests if a gene set is overrepresented in DE genes vs background
# Uses hypergeometric test (Fisher's exact test)

# Function to run ORA for a set of genes
run_ora <- function(genes, contrast_name, direction, output_path) {
  if (length(genes) == 0) {
    message(paste("No", direction, "genes for", contrast_name))
    return(NULL)
  }
  
  # Convert gene symbols to Entrez IDs
  # Required for most enrichment tools
  mapped_genes <- bitr(genes, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Mm.eg.db)
  sig_entrez <- mapped_genes$ENTREZID
  
  # Track unmapped genes
  unmapped_genes <- setdiff(genes, mapped_genes$SYMBOL)
  if (length(unmapped_genes) > 0) {
    message(paste("Unmapped", direction, "genes:", length(unmapped_genes)))
    write.csv(data.frame(symbol = unmapped_genes), 
             file.path(output_path, 
                      paste0(contrast_name, "_", direction, "_unmapped.csv")), 
             row.names = FALSE)
  }
  
  # ORA for GO Biological Process
  # Gene Ontology: hierarchical classification of gene functions
  ora_go <- enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",              # Biological Process (alt: MF, CC)
    pAdjustMethod = "BH",    # Benjamini-Hochberg FDR correction
    minGSSize = 10,          # Minimum genes per pathway
    pvalueCutoff = 0.05
  )
  
  # ORA for Reactome pathways
  # Reactome: curated database of biological pathways
  ora_reactome <- NULL
  try({
    ora_reactome <- enrichPathway(
      gene = sig_entrez,
      organism = "mouse",
      pAdjustMethod = "BH",
      minGSSize = 10,
      pvalueCutoff = 0.05,
      readable = TRUE         # Convert Entrez to symbols in output
    )
  }, silent = TRUE)
  
  # ORA for KEGG pathways
  ora_kegg <- NULL
  try({
    ora_kegg <- enricher(
      gene = sig_entrez,
      TERM2GENE = kegg_gmt,
      pAdjustMethod = "BH",
      minGSSize = 10,
      pvalueCutoff = 0.05
    )
  }, silent = TRUE)
  
  # Save results to CSV
  if (!is.null(ora_go) && nrow(as.data.frame(ora_go)) > 0) {
    write.csv(as.data.frame(ora_go), 
             file.path(output_path, 
                      paste0(contrast_name, "_", direction, "_ORA_GO.csv")), 
             row.names = FALSE)
  }
  if (!is.null(ora_reactome) && nrow(as.data.frame(ora_reactome)) > 0) {
    write.csv(as.data.frame(ora_reactome), 
             file.path(output_path, 
                      paste0(contrast_name, "_", direction, "_ORA_Reactome.csv")), 
             row.names = FALSE)
  }
  if (!is.null(ora_kegg) && nrow(as.data.frame(ora_kegg)) > 0) {
    write.csv(as.data.frame(ora_kegg), 
             file.path(output_path, 
                      paste0(contrast_name, "_", direction, "_ORA_KEGG.csv")), 
             row.names = FALSE)
  }
  
  # Create dotplots
  # Dotplot shows enriched terms with:
  # - x-axis: gene ratio (proportion of DE genes in pathway)
  # - Color: adjusted p-value
  # - Size: number of genes
  
  if (!is.null(ora_go) && nrow(as.data.frame(ora_go)) > 0) {
    p <- dotplot(ora_go, showCategory = 20, font.size = 7, 
                title = paste("GO BP:", contrast_name, direction))
    ggsave(file.path(output_path, 
                    paste0(contrast_name, "_", direction, "_ORA_GO_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(ora_reactome) && nrow(as.data.frame(ora_reactome)) > 0) {
    p <- dotplot(ora_reactome, showCategory = 20, font.size = 7, 
                title = paste("Reactome:", contrast_name, direction))
    ggsave(file.path(output_path, 
                    paste0(contrast_name, "_", direction, "_ORA_Reactome_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(ora_kegg) && nrow(as.data.frame(ora_kegg)) > 0) {
    p <- dotplot(ora_kegg, showCategory = 20, font.size = 7, 
                title = paste("KEGG:", contrast_name, direction))
    ggsave(file.path(output_path, 
                    paste0(contrast_name, "_", direction, "_ORA_KEGG_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  return(list(GO_BP = ora_go, Reactome = ora_reactome, KEGG = ora_kegg))
}

# Perform ORA for each contrast
ora_results <- list()
for (contrast_name in names(contrasts)) {
  cat("\n=== Running ORA for:", contrast_name, "===\n")
  
  # Extract significant genes
  sig_genes_df <- res_shrunken_sig_list[[contrast_name]]
  if (nrow(sig_genes_df) == 0) {
    message(paste("No significant genes for", contrast_name))
    next
  }
  
  # Split into up- and down-regulated genes
  # Analyze separately to identify direction-specific pathways
  up_genes <- sig_genes_df %>% filter(log2FoldChange > 0) %>% pull(symbol)
  down_genes <- sig_genes_df %>% filter(log2FoldChange < 0) %>% pull(symbol)
  
  # Run ORA for up-regulated genes
  ora_results[[contrast_name]][["up"]] <- 
    run_ora(up_genes, contrast_name, "up", output_path_ora)
  
  # Run ORA for down-regulated genes
  ora_results[[contrast_name]][["down"]] <- 
    run_ora(down_genes, contrast_name, "down", output_path_ora)
}

## GENE SET ENRICHMENT ANALYSIS (GSEA) ----

# GSEA uses all genes (not just significant ones)
# Tests if genes in a pathway are enriched at top/bottom of ranked list
# More powerful than ORA as it uses continuous rankings

# Function to run GSEA
run_gsea <- function(gene_list, contrast_name, output_path) {
  if (length(gene_list) == 0) {
    message(paste("No genes for GSEA:", contrast_name))
    return(NULL)
  }
  
  # GSEA for GO Biological Process
  gsea_go <- NULL
  try({
    gsea_go <- gseGO(
      geneList = gene_list,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      minGSSize = 10,
      pvalueCutoff = 0.05,
      eps = 0,               # Permutation-based p-value
      seed = 123
    )
  }, silent = FALSE)
  
  # GSEA for Reactome pathways
  gsea_reactome <- NULL
  try({
    gsea_reactome <- gsePathway(
      geneList = gene_list,
      organism = "mouse",
      pAdjustMethod = "BH",
      minGSSize = 10,
      pvalueCutoff = 0.05,
      eps = 0,
      seed = 123,
      verbose = FALSE
    )
  }, silent = FALSE)
  
  # GSEA for KEGG pathways
  gsea_kegg <- NULL
  try({
    gsea_kegg <- gseKEGG(
      geneList = gene_list,
      organism = "mmu",      # Mouse (change for other organisms)
      pAdjustMethod = "BH",
      minGSSize = 10,
      pvalueCutoff = 0.05,
      eps = 0,
      seed = 123,
      verbose = FALSE
    )
  }, silent = FALSE)
  
  # Save results
  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    write.csv(as.data.frame(gsea_go), 
             file.path(output_path, paste0(contrast_name, "_GSEA_GO.csv")), 
             row.names = FALSE)
  }
  if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
    write.csv(as.data.frame(gsea_reactome), 
             file.path(output_path, paste0(contrast_name, "_GSEA_Reactome.csv")), 
             row.names = FALSE)
  }
  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    write.csv(as.data.frame(gsea_kegg), 
             file.path(output_path, paste0(contrast_name, "_GSEA_KEGG.csv")), 
             row.names = FALSE)
  }
  
  # Create dotplots
  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    p <- dotplot(gsea_go, showCategory = 20, font.size = 7, 
                title = paste("GSEA GO:", contrast_name))
    ggsave(file.path(output_path, paste0(contrast_name, "_GSEA_GO_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
    p <- dotplot(gsea_reactome, showCategory = 20, font.size = 7, 
                title = paste("GSEA Reactome:", contrast_name))
    ggsave(file.path(output_path, paste0(contrast_name, "_GSEA_Reactome_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    p <- dotplot(gsea_kegg, showCategory = 20, font.size = 7, 
                title = paste("GSEA KEGG:", contrast_name))
    ggsave(file.path(output_path, paste0(contrast_name, "_GSEA_KEGG_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  return(list(GO_BP = gsea_go, Reactome = gsea_reactome, KEGG = gsea_kegg))
}

# Perform GSEA for each contrast
gsea_results <- list()
for (contrast_name in names(contrasts)) {
  cat("\n=== Running GSEA for:", contrast_name, "===\n")
  
  # Prepare ranked gene list
  # Ranking metric: sign(log2FC) * -log10(padj)
  # This ranks genes by both magnitude and significance
  gene_data <- res_shrunken_all_list[[contrast_name]] %>%
    filter(!is.na(padj)) %>%
    mutate(rank = sign(log2FoldChange) * -log10(padj))
  
  # Convert symbols to Entrez IDs
  input_symbols <- gene_data$symbol
  entrez_map <- bitr(input_symbols, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Mm.eg.db)
  
  # Save unmapped genes
  unmapped <- input_symbols[!input_symbols %in% entrez_map$SYMBOL]
  if (length(unmapped) > 0) {
    message(paste("Unmapped genes:", length(unmapped), 
                 "(", round(length(unmapped)/length(input_symbols)*100, 2), "%)"))
    write.csv(data.frame(symbol = unmapped), 
             file.path(output_path_gsea, paste0(contrast_name, "_unmapped.csv")), 
             row.names = FALSE)
  }
  
  # Filter to mapped genes
  gene_data <- gene_data %>%
    inner_join(entrez_map, by = c("symbol" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)  # Remove duplicate Entrez IDs
  
  # Create ranked gene list
  # Genes must be sorted in descending order of rank
  gene_list <- gene_data %>%
    arrange(desc(rank)) %>%
    pull(rank, name = ENTREZID)
  
  # Run GSEA
  gsea_results[[contrast_name]] <- run_gsea(gene_list, contrast_name, output_path_gsea)
}


################################################################################
# EXAMPLE USAGE
################################################################################

# To run the complete pipeline:
#
# 1. Prepare input files:
#    - data/gene_counts.tsv: Raw count matrix (genes x samples)
#    - data/sample_metadata.tsv: Sample information (sample, condition, ...)
#
# 2. Run QC and filtering (Script 1):
#    source("01_qc.R")
#
# 3. Run differential expression (Script 2):
#    source("02_de.R")
#
# 4. Run functional enrichment (Script 3):
#    source("03_enrichment.R")
#
# Output structure:
#   results/
#   ├── qc/          # QC plots and filtered data
#   ├── de/          # DE results and visualizations
#   ├── ora/         # Overrepresentation analysis results
#   └── gsea/        # GSEA results


################################################################################
# NOTES
################################################################################

# For different organisms:
# - Change org.Mm.eg.db to org.Hs.eg.db (human), org.Dm.eg.db (fly), etc.
# - Change msigdbr species to match
# - Change ReactomePA organism parameter
# - Change gseKEGG organism code (mmu = mouse, hsa = human, etc.)
#
# Common filtering thresholds:
# - Low expression: filterByExpr (adaptive) or rowSums(counts >= 10) >= 3
# - Significance: padj < 0.05 (FDR control)
# - Effect size: |log2FC| > 1 (2-fold change) or > 0.585 (1.5-fold)
#
# Alternative ranking metrics for GSEA:
# - Option 1: sign(log2FC) * -log10(padj)  # Recommended
# - Option 2: sign(log2FC) * -log10(pvalue)  # Without multiple testing correction
# - Option 3: log2FC  # Magnitude only (less common)
#
# Shrinkage methods:
# - apeglm: Adaptive prior (recommended, requires coef)
