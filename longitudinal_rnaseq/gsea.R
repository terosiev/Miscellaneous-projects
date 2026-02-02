# Gene Set Enrichment Analysis (GSEA)
# Author: Tero Sievänen
# Description: GSEA workflow using MSigDB gene sets for longitudinal RNA-seq data

# Load libraries ----
library(AnnotationHub)
library(clusterProfiler)
library(cowplot)
library(enrichplot)
library(tidyverse)
library(ensembldb)
library(msigdbr)
library(org.Hs.eg.db)
library(pathview)
library(grid)
library(ggplotify)

# Define paths ----
# NOTE: Update this path to match your directory structure
output_path_gsea <- "results/gsea/"



# Setup ----

# Get gene set collection from MSigDb
terms <- msigdbr(species = "Homo sapiens") %>% 
  dplyr::select(gs_name, gene_symbol)

# Define which contrasts to analyze
contrasts <- names(res_all_list)



# GSEA analysis ----

# Loop through contrasts
gsea_results <- list()

for (contrast in contrasts) {
  
  message(paste("Running GSEA for contrast:", contrast))
  
  # Get all genes for this contrast (not just significant)
  contrast_data <- res_all_list[[contrast]]
  
  # Sort by logFC
  contrast_data <- contrast_data[order(contrast_data$logFC, decreasing = TRUE), ]
  
  # Extract the log2 fold change
  gene_list <- contrast_data$logFC
  
  # Name the vector using gene symbols
  names(gene_list) <- contrast_data$symbol
  
  # Remove NA symbols
  gene_list <- gene_list[!is.na(names(gene_list))]
  
  # Remove duplicates (keep first occurrence - highest logFC due to sorting)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  
  # Check we have enough genes
  if (length(gene_list) < 100) {
    message(paste("  Skipping", contrast, "- too few genes (n =", length(gene_list), ")"))
    next
  }
  
  # Inspect
  message(paste("  Gene list length:", length(gene_list)))
  message(paste("  Top gene:", names(gene_list)[1], "=", gene_list[1]))
  message(paste("  Bottom gene:", names(gene_list)[length(gene_list)], "=", gene_list[length(gene_list)]))
  
  # Run GSEA
  gse <- GSEA(geneList = gene_list, 
              TERM2GENE = terms, 
              pvalueCutoff = 0.05, 
              pAdjustMethod = "BH", 
              nPermSimple = 10000, 
              minGSSize = 10)
  
  # Store results
  gsea_results[[contrast]] <- gse
  
  # Save results
  gse.res <- data.frame(gse@result)
  write.csv(gse.res, file = file.path(output_path_gsea, paste0(contrast, "_gsea.csv")), 
            quote = FALSE, row.names = FALSE)
  
  message(paste("  Found", nrow(gse.res), "enriched gene sets"))
}



# Visualization ----

message("Creating GSEA visualizations...")

for (contrast in names(gsea_results)) {
  
  gse <- gsea_results[[contrast]]
  
  # Skip if no results
  if (is.null(gse) || nrow(gse) == 0) {
    message(paste("  Skipping visualization for", contrast, "- no enriched gene sets"))
    next
  }
  
  # Get gene list for this contrast
  contrast_data <- res_all_list[[contrast]]
  contrast_data <- contrast_data[order(contrast_data$logFC, decreasing = TRUE), ]
  gene_list <- contrast_data$logFC
  names(gene_list) <- contrast_data$symbol
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- gene_list[!duplicated(names(gene_list))]
  
  # Plot dot plot for enriched gene sets
  pdf(file.path(output_path_gsea, paste0(contrast, "_gsea_dotplot.pdf")), width = 9, height = 12.5)
  print(dotplot(gse, showCategory=20, split=".sign", font.size = 8) + facet_grid(.~.sign))
  dev.off()
  
  # Plot top 10 enrichment plots
  n_plots <- min(10, nrow(gse))
  
  if (n_plots > 0) {
    pdf_list <- list()
    
    # Generate individual plots and store them in the list
    for (i in 1:n_plots) {
      pdf_list[[i]] <- as.grob(gseaplot2(gse, title = gse$Description[i], geneSetID = i))
    }
    
    # Create the grid plot
    grid_plot <- plot_grid(plotlist = pdf_list, ncol = 2)
    
    # Save the grid plot (5 rows x 2 columns = reasonable size)
    ggsave(file.path(output_path_gsea, paste0(contrast, "_gsea_grid_plot.pdf")), 
           grid_plot, width = 16, height = 20)
  }
  
  # Ridgeplot (only if we have results)
  if (nrow(gse) > 0) {
    pdf(file.path(output_path_gsea, paste0(contrast, "_gsea_ridge.pdf")), width = 16, height = 13)
    print(ridgeplot(gse, showCategory = min(20, nrow(gse)), label_format = 50) + 
            labs(x = "Enrichment distribution") + 
            ggtitle(paste("GSEA Enrichment distributions -", contrast)))
    dev.off()
  }
  
  # Extract top positively and negatively enriched gene sets
  if (sum(gse$enrichmentScore > 0) > 0) {
    topUpper <- gse$Description[gse$enrichmentScore > 0][1:min(5, sum(gse$enrichmentScore > 0))]
    
    # cnetplot up
    pdf(file.path(output_path_gsea, paste0(contrast, "_gsea_cnet_up.pdf")), width = 15, height = 8)
    print(cnetplot(gse, categorySize="geneNum", foldChange=gene_list, 
                   showCategory = topUpper, layout = "fr", colorEdge = TRUE, 
                   cex.params = list(category_label = 0.7, gene_label = 0.5)))
    dev.off()
  }
  
  if (sum(gse$enrichmentScore < 0) > 0) {
    topLower <- gse$Description[gse$enrichmentScore < 0][1:min(5, sum(gse$enrichmentScore < 0))]
    
    # cnetplot down
    pdf(file.path(output_path_gsea, paste0(contrast, "_gsea_cnet_low.pdf")), width = 15, height = 8)
    print(cnetplot(gse, categorySize="geneNum", foldChange=gene_list, 
                   showCategory = topLower, layout = "fr", colorEdge = TRUE, 
                   cex.params = list(category_label = 0.7, gene_label = 0.5)))
    dev.off()
  }
  
  message(paste("  Visualization complete for", contrast))
}

message("GSEA analysis complete!")

# Print session info
sessionInfo()
