################################################################################
# Script 1: Quality Control and Data Filtering
#
# This script performs:
# - Data import and inspection
# - Library size and count distribution visualization
# - PCA analysis (before and after filtering)
# - Low expression filtering
# - Sample correlation heatmaps
#
# Outputs: results/qc/
################################################################################

# Load configuration
source("scripts/00_config.R")

cat("\n=== Starting QC and Filtering ===\n")

## DATA IMPORT ----

cat("Importing data...\n")

# Import raw count matrix
counts.df <- data.frame(fread(COUNTS_FILE, sep = "\t"))

# Import sample metadata
meta <- data.frame(fread(METADATA_FILE, sep = "\t"))

# Store gene annotation
gene_annotation <- counts.df %>% 
  dplyr::select(gene_id)
gene_annotation$symbol <- counts.df$gene_name

# Prepare count matrix
rownames(counts.df) <- counts.df$gene_id
counts.df <- counts.df %>% 
  dplyr::select(-gene_id, -gene_name)
counts.df <- round(counts.df)

# Prepare metadata
rownames(meta) <- meta$sample

cat(sprintf("Loaded %d genes and %d samples\n", nrow(counts.df), ncol(counts.df)))

## DATA INSPECTION ----

cat("Calculating library sizes...\n")

# Library sizes
library_sizes <- data.frame(
  Sample = colnames(counts.df),
  TotalCounts = colSums(counts.df),
  Condition = meta$condition[match(colnames(counts.df), rownames(meta))]
)

mean_counts <- mean(library_sizes$TotalCounts)
median_counts <- median(library_sizes$TotalCounts)

# Plot library sizes
p <- ggplot(library_sizes, aes(x = Sample, y = TotalCounts, fill = Condition)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean_counts, linetype = "solid", color = "black", linewidth = 1) +
  geom_hline(yintercept = median_counts, linetype = "dashed", color = "gray40", linewidth = 1) +
  scale_fill_manual(values = ANN_COLORS$condition) +
  ggtitle("Library Sizes") +
  xlab("Sample") +
  ylab("Total Counts") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(QC_DIR, "library_sizes.pdf"), plot = p, width = 8, height = 6)

# Plot count distribution
count_dist <- counts.df %>%
  mutate(Gene = rownames(counts.df)) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Counts") %>%
  mutate(LogCounts = log2(Counts + 1)) %>%
  left_join(dplyr::select(meta, sample, condition), by = c("Sample" = "sample"))

p <- ggplot(count_dist, aes(x = Sample, y = LogCounts, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = ANN_COLORS$condition) +
  ggtitle("Log2 Count Distribution") +
  xlab("Sample") +
  ylab("Log2(Counts + 1)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(QC_DIR, "count_distribution.pdf"), plot = p, width = 8, height = 6)

## CREATE DESEQ2 OBJECT ----

cat("Creating DESeq2 object...\n")

colData <- meta %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "sample") %>%
  dplyr::select(condition)

dds <- DESeqDataSetFromMatrix(
  countData = counts.df, 
  colData = colData, 
  design = ~ condition
)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Plot dispersion
pdf(file.path(QC_DIR, "dispersion.pdf"))
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

## PCA FUNCTION ----

run_PCA <- function(dds, design_matrix, pdf_filename, path) {
  rld <- vst(dds, blind = TRUE)
  pca <- prcomp(t(assay(rld)))
  var_explained <- summary(pca)$importance[2, ] * 100
  
  df <- data.frame(cbind(design_matrix, pca$x))
  
  p <- ggplot(df, aes(x = PC1, y = PC2, color = condition)) + 
    geom_point(size = 7, alpha = 0.8) +
    geom_point(size = 7, shape = 1, color = "black") +
    geom_text(aes(label = rownames(design_matrix)), 
             nudge_x = 0.8, nudge_y = 0.8, 
             check_overlap = TRUE, color = "black") +
    theme_classic() +
    ggtitle("PCA Plot") +
    xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    scale_color_manual(values = ANN_COLORS$condition)
  
  ggsave(file.path(path, paste0(pdf_filename, ".pdf")), plot = p, width = 8, height = 6, dpi = 300)
  
  return(p)
}

## PCA PRE-FILTERING ----

cat("Running PCA (pre-filtering)...\n")
run_PCA(dds, meta, "PCA_prefiltering", QC_DIR)

## FILTERING ----

cat("Filtering low-expression genes...\n")

# Use edgeR's filterByExpr for adaptive filtering
keep <- filterByExpr(dds, design = model.matrix(~ condition, data = colData))
dds_filtered <- dds[keep, ]

cat(sprintf("Retained %d / %d genes (%.1f%%)\n", 
            sum(keep), length(keep), sum(keep)/length(keep)*100))

# Update gene annotation
gene_annotation_filtered <- gene_annotation[keep, ]

## PCA POST-FILTERING ----

cat("Running PCA (post-filtering)...\n")
run_PCA(dds_filtered, meta, "PCA_postfiltering", QC_DIR)

## SAMPLE CORRELATION ----

cat("Generating correlation heatmap...\n")

rld <- vst(dds_filtered, blind = TRUE)
sample_cor <- cor(assay(rld), method = "pearson")

annotation_col <- data.frame(
  condition = meta$condition,
  row.names = rownames(meta)
)

pdf(file.path(QC_DIR, "sample_correlation.pdf"), width = 10, height = 10)
pheatmap(sample_cor,
         annotation_col = annotation_col,
         annotation_colors = ANN_COLORS,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Sample Correlation (Pearson)")
dev.off()

## SAVE FILTERED DATA ----

cat("Saving filtered data...\n")

save(dds_filtered, gene_annotation_filtered, meta, 
     file = file.path(QC_DIR, "filtered_data.RData"))

cat("\n=== QC and Filtering Complete ===\n")
cat(sprintf("Outputs saved to: %s\n", QC_DIR))
