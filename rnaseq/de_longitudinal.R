# Longitudinal differential expression analysis
# Author: Tero Sievänen
# Description: DE analysis using limma-voom with duplicate correlation for longitudinal RNA-seq data

# Load libraries ----
library(DESeq2)
library(DOSE)
library(edgeR)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(tximport)
library(limma)
library(EnhancedVolcano)

# Define paths ----
# NOTE: Update this path to match your directory structure
output_path_deg <- "results/differential_expression/"

# Differential expression analysis ----

# Define variables (not mandatory)
condition <- data$samples$condition
time <- data$samples$time
batch <- data$samples$batch
group <- factor(paste(condition,time,sep="_"))

# Create the model matrix
design <- model.matrix(~0+ group + batch)

# Inspect
design

# Remove heteroscedasticity - first voom
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)

# Estimate the correlation between measurements made on the same subject
corfit <- duplicateCorrelation(v,design,block=data$samples$subject)

# Inspect
corfit$consensus

# Re-run voom with correlation 
v <- voom(data, design, plot=TRUE, block=data$samples$subject, correlation=corfit$consensus)

# This inter-subject correlation is input into the linear model fit:
fit <- lmFit(v, design, block=data$samples$subject, correlation=corfit$consensus)

# Create contrasts
# NOTE: Update these contrasts to match your experimental design
cm <- makeContrasts(
  Cond1_T2T1 = groupCONDITION1_time2 - groupCONDITION1_time1,
  Cond1_T3T1 = groupCONDITION1_time3 - groupCONDITION1_time1,
  Cond1_T3T2 = groupCONDITION1_time3 - groupCONDITION1_time2,
  Cond2_T2T1 = groupCONDITION2_time2 - groupCONDITION2_time1,
  Cond2_T3T1 = groupCONDITION2_time3 - groupCONDITION2_time1,
  Cond2_T3T2 = groupCONDITION2_time3 - groupCONDITION2_time2,
  Cond1_Cond2_T1 = groupCONDITION1_time1 - groupCONDITION2_time1,
  Cond1_Cond2_T2 = groupCONDITION1_time2 - groupCONDITION2_time2,
  Cond1_Cond2_T3 = groupCONDITION1_time3 - groupCONDITION2_time3,
  MainEffect = (groupCONDITION1_time3 + groupCONDITION1_time2 + groupCONDITION1_time1)/3 - (groupCONDITION2_time3 + groupCONDITION2_time2 + groupCONDITION2_time1)/3,
  Interaction_T3T1 = (groupCONDITION1_time3 - groupCONDITION1_time1) - (groupCONDITION2_time3 - groupCONDITION2_time1),
  Interaction_T2T1 = (groupCONDITION1_time2 - groupCONDITION1_time1) - (groupCONDITION2_time2 - groupCONDITION2_time1),
  Interaction_T3T2 = (groupCONDITION1_time3 - groupCONDITION1_time2) - (groupCONDITION2_time3 - groupCONDITION2_time2),
  levels=design)

# Fit contrasts
fit2 <- contrasts.fit(fit, cm)

# Run DE analysis with limma
efit <- eBayes(fit2)

# Plot
plotSA(efit, main="Final model: Mean-variance trend")

# Inspect
summary(decideTests(efit))

# Store all results
res_all <- topTable(efit, n = Inf)
colnames(res_all)[1] <- "symbol"
res_all$gene <- rownames(res_all)

# Define the contrasts
contrasts <- c("Cond1_T2T1", "Cond1_T3T1", "Cond1_T3T2", 
               "Cond2_T2T1", "Cond2_T3T1", "Cond2_T3T2",
               "Cond1_Cond2_T1", "Cond1_Cond2_T2", "Cond1_Cond2_T3",
               "MainEffect", "Interaction_T2T1", "Interaction_T3T1", "Interaction_T3T2")


# Initialize a list to store the results
res_sig_list <- list()
res_all_list <- list()

# Loop through the contrasts and print DEG lists
for (contrast in contrasts) {
  # Extract topTable results for the current contrast
  res_sig <- topTable(efit, coef = contrast, n = Inf, sort = "p", p = 0.05)[,-1]
  res_all_c <- topTable(efit, coef = contrast, n = Inf, sort = "p")
  
  # Add gene labels to results table
  res_sig$gene <- rownames(res_sig)
  res_sig$symbol <- grch38annot$symbol[match(rownames(res_sig), grch38annot$gene)]
  res_all_c$gene <- rownames(res_all_c)
  res_all_c$symbol <- grch38annot$symbol[match(rownames(res_all_c), grch38annot$gene)]
  
  # Store the results in the list
  res_sig_list[[contrast]] <- res_sig
  res_all_list[[contrast]] <- res_all_c
  
  # Save significant DEGs
  write.csv(as.data.frame(res_sig_list[[contrast]]), 
            file = file.path(output_path_deg, paste0(contrast, "_DEG_list.csv")), 
            row.names = TRUE)
  
  # Save full results for all genes
  write.csv(as.data.frame(res_all_list[[contrast]]), 
            file = file.path(output_path_deg, paste0(contrast, "_all_genes.csv")), 
            row.names = TRUE)
}

# Data visualizations ----

# Draw volcano plots of each contrast
for (contrast in contrasts) {
  # Extract data for the current contrast
  contrast_data <- res_all_list[[contrast]]
  
  # Set plot width and height
  options(repr.plot.width=20, repr.plot.height=20)
  
  # Check if there are any significant genes
  n_sig <- sum(contrast_data$adj.P.Val < 0.05 & abs(contrast_data$logFC) > 1, na.rm = TRUE)
  
  # Adjust parameters based on whether there are significant hits
  if (n_sig > 0) {
    # Normal volcano plot with labeling
    volcano_plot <- EnhancedVolcano(contrast_data, x = "logFC",
                                    y = "adj.P.Val", 
                                    ylim = c(0, max(-log10(contrast_data[["adj.P.Val"]]), na.rm = TRUE) + 0.5), 
                                    lab = contrast_data$symbol,
                                    labSize = 3,
                                    colAlpha = 0.8, 
                                    pCutoff = 0.05,
                                    FCcutoff = 1,
                                    col = c('grey30', '#66C2A5', '#8DA0CB', '#FC8D62'),
                                    title = paste("DE genes", contrast),
                                    subtitle = paste(n_sig, "significant genes"),
                                    legendPosition = "bottom",
                                    legendLabels = c("NS", "Log2 FC", "P-value", "Sig"),
                                    drawConnectors = FALSE) +
      theme_classic()
  } else {
    # No significant hits - simpler plot without labels
    volcano_plot <- EnhancedVolcano(contrast_data, x = "logFC",
                                    y = "adj.P.Val", 
                                    ylim = c(0, max(-log10(contrast_data[["adj.P.Val"]]), na.rm = TRUE) + 0.5), 
                                    lab = NA,  # No labels
                                    colAlpha = 0.8, 
                                    pCutoff = 0.05,
                                    FCcutoff = 1,
                                    col = c('grey30', '#66C2A5', '#8DA0CB', '#FC8D62'),
                                    title = paste("DE genes", contrast),
                                    subtitle = "No significant genes (adj.P < 0.05 & |logFC| > 1)",
                                    legendPosition = "bottom",
                                    legendLabels = c("NS", "Log2 FC", "P-value", "P & FC"),
                                    drawConnectors = FALSE) +
      theme_classic()
  }
  
  # Save the volcano plot as a PDF with DPI 300
  ggsave(filename = file.path(output_path_deg, paste0(contrast, "_volcano_plot.pdf")), 
         plot = volcano_plot, dpi = 300)
  
}

# Extract normalized counts
normalized_counts <- v$E %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  left_join(grch38annot, by = c("gene" = "gene")) %>%
  dplyr::select(gene, symbol, everything())

# Draw heatmaps of each contrast
for (contrast in contrasts) {
  # Extract normalized expression for significant genes
  norm.sig <- normalized_counts %>% 
    dplyr::filter(gene %in% rownames(res_sig_list[[contrast]]))
  
  # Only plot if there are at least 2 significant genes
  if (nrow(norm.sig) < 2) {
    message(paste("Skipping heatmap for", contrast, "- fewer than 2 significant genes (n =", nrow(norm.sig), ")"))
    next
  }
  
  # Determine the number of genes to display (up to 100 at maximum)
  num_genes <- min(nrow(norm.sig), 100)
  
  # Subset the top genes
  top_genes <- norm.sig[1:num_genes, ]
  rownames(top_genes) <- top_genes$symbol
  
  # Calculate z-scores for each gene across samples
  top_genes_zscore <- data.frame(t(apply(top_genes[, 3:ncol(top_genes)], 1, scale)))
  colnames(top_genes_zscore) <- colnames(counts)
  
  # Choose annotations
  anno <- targets[, -c(3,4)]
  
  # Set a color palette
  heat_colors <- rev(brewer.pal(7, "RdBu"))
  
  # Define annotation colors
  ann_colors = list(
    time = c(time1 ="#66C2A5", 
             time2="#FC8D62",
             time3="#8DA0CB"),
    condition = c(CONDITION1 = "#D53E4F", CONDITION2 = "#3288BD")
  )
  
  # Plot DEG heatmap
  pheatmap(top_genes_zscore,
           color = heat_colors,
           cluster_rows = T,
           cluster_cols = F,
           show_rownames = T,
           annotation = anno, 
           annotation_colors = ann_colors,
           fontsize = 5,
           fontsize_row = 5,
           angle_col = 45,
           filename = file.path(output_path_deg, paste0(contrast, "_hmap", ".pdf")))
  
  message(paste("Heatmap created for", contrast, "with", num_genes, "genes"))
}
