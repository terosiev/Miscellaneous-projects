# Over-representation analysis (ORA)
# Author: Tero Sievänen
# Description: Gene ontology enrichment analysis for differentially expressed genes

# Load libraries ----
library(AnnotationHub)
library(clusterProfiler)
library(tidyverse)
library(ensembldb)
library(msigdbr)
library(org.Hs.eg.db)

# Define paths ----
# NOTE: Update this path to match your directory structure
output_path_ora <- "results/ora/"

# Define functions ----

# Gene ontology enrichment analysis
run_GO_enrichment <- function(genes, universe, file_name) {
  ego <- enrichGO(gene = genes, 
                  universe = universe,
                  keyType = "ENSEMBL",
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.01,
                  minGSSize = 20,
                  readable = TRUE)
  
  # Write ego to file
  write.csv(as.data.frame(ego), file = file.path(output_path_ora, file_name))
  
  # Return ego
  return(ego)
}

# Draw dotplots and cnetplots
enhanced_plot <- function(ego, file_name, output_path, OE_foldchanges, showCategory = 10) {
  if (!is.null(ego) && nrow(ego) > 0) {
    # Dotplot
    dotplot(ego, showCategory = showCategory, font.size = 7)
    ggsave(file.path(output_path, paste0("dotplot_", file_name, ".pdf")), units = "in", height = 12.5, width = 9)
    
    # Cnetplot
    cnetplot(ego, cex.params = list(category_label = 0.7, gene_label = 0.5), 
             color.params = list(foldChange = OE_foldchanges), 
             categorySize = "pvalue", showCategory = 5, colorEdge = T)
    ggsave(file.path(output_path, paste0("cnetplot_", file_name, ".pdf")), height = 10, width = 15)
  } else {
    message(paste("No enrichment results to plot for", file_name))
  }
}



# Setup ----

# Create a cache
ah <- AnnotationHub()

# Query AnnotationHub
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))

# Extract annotations of interest
# NOTE: Update the ID based on your genome version
human_ens <- human_ens[["AH116291"]]



# ORA analysis ----

# NOTE: This section should be run AFTER the DE analysis in de.R
# Required objects: data (filtered DGEList), res_sig_list, grch38annot

# Create background dataset using ALL genes that passed filtering
all_genes <- rownames(data) %>% sub("\\..*", "", .)

# Create a gene-level data frame with all filtered genes for annotation
annotations_ahb <- genes(human_ens, return.type = "data.frame") %>%
  dplyr::select(gene_id_version, gene_name, entrezid, gene_biotype) %>%
  dplyr::filter(gene_id_version %in% rownames(data))

# Choose only the first match for entrezid
annotations_ahb$entrezid <- map(annotations_ahb$entrezid, 1) %>% unlist()

# Remove duplicated gene names (keep first occurrence)
annotations_ahb <- annotations_ahb[!duplicated(annotations_ahb$gene_name), ]

# Remove NA entries in Entrez column
annotations_ahb <- annotations_ahb[!is.na(annotations_ahb$entrezid), ]

# Define which contrasts to analyze
contrasts <- names(res_sig_list)

# Loop through each contrast and perform ORA
ego_list <- list()

for (contrast in contrasts) {
  
  message(paste("Processing ORA for contrast:", contrast))
  
  # Get DEGs for this contrast
  DEG_df <- res_sig_list[[contrast]]
  
  # Skip if no significant genes
  if (nrow(DEG_df) == 0) {
    message(paste("  Skipping", contrast, "- no significant genes"))
    next
  }
  
  # Remove genes with padj = NA
  DEG_df_noNAs <- dplyr::filter(DEG_df, !is.na(adj.P.Val))
  
  # Merge the AnnotationHub dataframe with the results 
  res_ids <- left_join(DEG_df_noNAs, annotations_ahb, by = c("gene" = "gene_id_version"))
  
  # Remove rows with NA entrezid (couldn't be annotated)
  res_ids <- res_ids %>% dplyr::filter(!is.na(entrezid))
  
  # Extract significant gene IDs (Ensembl without version)
  sig_genes <- res_ids$gene %>% sub("\\..*", "", .)
  
  # Skip if too few genes after filtering
  if (length(sig_genes) < 10) {
    message(paste("  Skipping", contrast, "- fewer than 10 annotated genes"))
    next
  }
  
  # Separate upregulated and downregulated genes
  upregulated_genes <- res_ids %>%
    dplyr::filter(logFC > 0) %>%
    pull(gene) %>%
    sub("\\..*", "", .)
  
  downregulated_genes <- res_ids %>%
    dplyr::filter(logFC < 0) %>%
    pull(gene) %>%
    sub("\\..*", "", .)
  
  # Create fold change vector for plotting
  OE_foldchanges <- res_ids$logFC
  names(OE_foldchanges) <- res_ids$gene %>% sub("\\..*", "", .)
  
  # Run GO enrichment for all DEGs
  message(paste("  Running GO enrichment for all DEGs (n =", length(sig_genes), ")"))
  ego <- run_GO_enrichment(sig_genes, all_genes, 
                           paste0(contrast, "_all_DEGs.csv"))
  ego_list[[paste0(contrast, "_all")]] <- ego
  
  # Run GO enrichment for upregulated genes (if enough genes)
  if (length(upregulated_genes) >= 10) {
    message(paste("  Running GO enrichment for upregulated genes (n =", length(upregulated_genes), ")"))
    ego_up <- run_GO_enrichment(upregulated_genes, all_genes, 
                                paste0(contrast, "_upregulated.csv"))
    ego_list[[paste0(contrast, "_up")]] <- ego_up
  }
  
  # Run GO enrichment for downregulated genes (if enough genes)
  if (length(downregulated_genes) >= 10) {
    message(paste("  Running GO enrichment for downregulated genes (n =", length(downregulated_genes), ")"))
    ego_down <- run_GO_enrichment(downregulated_genes, all_genes, 
                                  paste0(contrast, "_downregulated.csv"))
    ego_list[[paste0(contrast, "_down")]] <- ego_down
  }
  
  # Store fold changes for visualization
  assign(paste0("OE_foldchanges_", contrast), OE_foldchanges, envir = .GlobalEnv)
}



# Visualization ----

message("Creating visualization plots...")

# Run enhanced plots for each contrast
for (ego_name in names(ego_list)) {
  
  ego <- ego_list[[ego_name]]
  
  # Skip if NULL or empty
  if (is.null(ego) || nrow(ego) == 0) {
    next
  }
  
  # Calculate pairwise term similarity for plotting
  ego <- enrichplot::pairwise_termsim(ego)
  
  # Extract contrast name to get the correct fold changes
  contrast_base <- sub("_(all|up|down)$", "", ego_name)
  OE_foldchanges <- get(paste0("OE_foldchanges_", contrast_base))
  
  # Create plots
  enhanced_plot(ego, ego_name, output_path_ora, OE_foldchanges)
  
  message(paste("  Plots created for", ego_name))
}

message("ORA analysis complete!")
