################################################################################
# Script 2: Differential Expression Analysis
#
# This script performs:
# - DESeq2 differential expression analysis
# - Log fold change shrinkage
# - Volcano plots and MA plots
# - Heatmaps of top DEGs
# - Export of DE results
#
# Outputs: results/de/
################################################################################

# Load configuration
source("scripts/00_config.R")

cat("\n=== Starting Differential Expression Analysis ===\n")

## LOAD FILTERED DATA ----

cat("Loading filtered data...\n")
load(file.path(QC_DIR, "filtered_data.RData"))

## RUN DESEQ2 ----

cat("Running DESeq2...\n")
dds_filtered <- DESeq(dds_filtered)

cat("Available coefficients:\n")
print(resultsNames(dds_filtered))

## DEFINE CONTRASTS ----

# Modify these contrasts to match your experimental design
contrasts <- list(
  Treatment1_vs_Control = c("condition", "Treatment1", "Control"),
  Treatment2_vs_Control = c("condition", "Treatment2", "Control")
)

## PERFORM DE ANALYSIS ----

res_list <- list()
de_summary <- data.frame(
  Contrast = character(),
  Total_Significant = integer(),
  Up_Regulated = integer(),
  Down_Regulated = integer(),
  stringsAsFactors = FALSE
)

for (contrast_name in names(contrasts)) {
  cat(sprintf("\n=== Analyzing: %s ===\n", contrast_name))
  
  # Extract results
  res <- results(dds_filtered, contrast = contrasts[[contrast_name]], 
                alpha = PADJ_THRESHOLD)
  
  # Shrink log2 fold changes (using ashr - works with all contrasts)
  res_shrunken <- lfcShrink(dds_filtered, 
                            contrast = contrasts[[contrast_name]], 
                            type = "ashr", 
                            res = res)
  
  # Add gene annotation
  res_df <- as.data.frame(res_shrunken) %>%
    mutate(gene_id = rownames(res_shrunken)) %>%
    left_join(gene_annotation_filtered, by = c("gene_id" = "gene_id")) %>%
    arrange(padj)
  
  # Filter significant genes
  res_sig <- res_df %>%
    filter(padj < PADJ_THRESHOLD & abs(log2FoldChange) > LFC_THRESHOLD)
  
  # Summary statistics
  total_sig <- nrow(res_sig)
  up_reg <- sum(res_sig$log2FoldChange > 0, na.rm = TRUE)
  down_reg <- sum(res_sig$log2FoldChange < 0, na.rm = TRUE)
  
  cat(sprintf("Significant genes: %d (Up: %d, Down: %d)\n", 
              total_sig, up_reg, down_reg))
  
  # Store results
  res_list[[contrast_name]] <- list(
    all = res_df,
    sig = res_sig,
    results_obj = res_shrunken
  )
  
  # Update summary table
  de_summary <- rbind(de_summary, data.frame(
    Contrast = contrast_name,
    Total_Significant = total_sig,
    Up_Regulated = up_reg,
    Down_Regulated = down_reg
  ))
  
  # Save CSV files
  write.csv(res_df, 
           file.path(DE_DIR, paste0(contrast_name, "_all_genes.csv")), 
           row.names = FALSE)
  
  if (nrow(res_sig) > 0) {
    write.csv(res_sig, 
             file.path(DE_DIR, paste0(contrast_name, "_sig_genes.csv")), 
             row.names = FALSE)
  }
}

# Save summary
write.csv(de_summary, file.path(DE_DIR, "de_summary.csv"), row.names = FALSE)
print(de_summary)

## VOLCANO PLOTS ----

cat("\nGenerating volcano plots...\n")

for (contrast_name in names(contrasts)) {
  res_df <- res_list[[contrast_name]]$all
  
  # Add significance category
  res_df <- res_df %>%
    mutate(
      sig_category = case_when(
        padj < PADJ_THRESHOLD & log2FoldChange > LFC_THRESHOLD ~ "Up",
        padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Count categories
  sig_counts <- table(res_df$sig_category)
  
  # Create volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig_category)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "gray70")) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), 
              linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(PADJ_THRESHOLD), 
              linetype = "dashed", color = "black") +
    labs(title = gsub("_", " ", contrast_name),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value",
         color = "Regulation") +
    annotate("text", x = Inf, y = Inf, 
            label = sprintf("Up: %d\nDown: %d", 
                          sig_counts["Up"], sig_counts["Down"]),
            hjust = 1.1, vjust = 1.1, size = 4) +
    theme_classic() +
    theme(legend.position = "top")
  
  ggsave(file.path(DE_DIR, paste0(contrast_name, "_volcano.pdf")), 
        plot = p, width = 8, height = 6, dpi = 300)
}

## MA PLOTS ----

cat("Generating MA plots...\n")

for (contrast_name in names(contrasts)) {
  pdf(file.path(DE_DIR, paste0(contrast_name, "_MA_plot.pdf")), width = 8, height = 6)
  plotMA(res_list[[contrast_name]]$results_obj, 
         ylim = c(-5, 5),
         main = gsub("_", " ", contrast_name))
  dev.off()
}

## HEATMAP OF TOP DEGS ----

cat("Generating heatmaps...\n")

for (contrast_name in names(contrasts)) {
  res_sig <- res_list[[contrast_name]]$sig
  
  if (nrow(res_sig) < 2) {
    cat(sprintf("Skipping heatmap for %s (< 2 significant genes)\n", contrast_name))
    next
  }
  
  # Get top 50 significant genes (or all if fewer)
  top_genes <- head(res_sig$gene_id, 50)
  
  # Extract normalized counts
  vsd <- vst(dds_filtered, blind = FALSE)
  mat <- assay(vsd)[top_genes, ]
  
  # Scale rows (z-score)
  mat_scaled <- t(scale(t(mat)))
  
  # Create annotation
  annotation_col <- data.frame(
    Condition = colData(dds_filtered)$condition,
    row.names = colnames(mat_scaled)
  )
  
  # Row labels (gene symbols)
  gene_labels <- gene_annotation_filtered %>%
    filter(gene_id %in% top_genes) %>%
    arrange(match(gene_id, top_genes))
  rownames(mat_scaled) <- gene_labels$symbol
  
  # Plot heatmap
  pdf(file.path(DE_DIR, paste0(contrast_name, "_heatmap_top50.pdf")), 
      width = 10, height = 12)
  pheatmap(mat_scaled,
           annotation_col = annotation_col,
           annotation_colors = ANN_COLORS,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
           breaks = seq(-2, 2, length.out = 101),
           main = paste("Top 50 DEGs:", gsub("_", " ", contrast_name)),
           fontsize_row = 8,
           show_colnames = TRUE,
           cluster_cols = TRUE,
           cluster_rows = TRUE)
  dev.off()
}

## SAVE WORKSPACE ----

save(res_list, de_summary, dds_filtered, gene_annotation_filtered,
     file = file.path(DE_DIR, "de_results.RData"))

cat("\n=== Differential Expression Analysis Complete ===\n")
cat(sprintf("Outputs saved to: %s\n", DE_DIR))
