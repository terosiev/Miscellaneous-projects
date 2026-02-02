################################################################################
# Script 3: Functional Enrichment Analysis
#
# This script performs:
# - Over-Representation Analysis (ORA) for up/down-regulated genes
# - Gene Set Enrichment Analysis (GSEA) using ranked gene lists
# - Pathway analysis using GO, Reactome, and KEGG databases
#
# Outputs: results/ora/ and results/gsea/
################################################################################

# Load configuration
source("scripts/00_config.R")

cat("\n=== Starting Functional Enrichment Analysis ===\n")

## LOAD DE RESULTS ----

cat("Loading DE results...\n")
load(file.path(DE_DIR, "de_results.RData"))

# Get organism database
orgdb <- get(ORGDB)

## OVER-REPRESENTATION ANALYSIS (ORA) ----

cat("\n=== Running Over-Representation Analysis ===\n")

run_ora <- function(gene_symbols, contrast_name, direction, output_path) {
  
  if (length(gene_symbols) == 0) {
    message(sprintf("No %s genes for %s", direction, contrast_name))
    return(NULL)
  }
  
  cat(sprintf("\nAnalyzing %d %s genes for %s\n", 
              length(gene_symbols), direction, contrast_name))
  
  # Convert symbols to Entrez IDs
  entrez_map <- bitr(gene_symbols, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = orgdb)
  
  gene_entrez <- entrez_map$ENTREZID
  
  if (length(gene_entrez) == 0) {
    message("No genes mapped to Entrez IDs")
    return(NULL)
  }
  
  # GO Biological Process
  ora_go <- NULL
  try({
    ora_go <- enrichGO(
      gene = gene_entrez,
      OrgDb = orgdb,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = ENRICHMENT_PADJ,
      qvalueCutoff = ENRICHMENT_PADJ,
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE,
      readable = TRUE
    )
  }, silent = FALSE)
  
  # Reactome
  ora_reactome <- NULL
  try({
    ora_reactome <- enrichPathway(
      gene = gene_entrez,
      organism = ORGANISM,
      pAdjustMethod = "BH",
      pvalueCutoff = ENRICHMENT_PADJ,
      qvalueCutoff = ENRICHMENT_PADJ,
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE,
      readable = TRUE
    )
  }, silent = FALSE)
  
  # KEGG
  ora_kegg <- NULL
  try({
    ora_kegg <- enrichKEGG(
      gene = gene_entrez,
      organism = KEGG_ORGANISM,
      pAdjustMethod = "BH",
      pvalueCutoff = ENRICHMENT_PADJ,
      qvalueCutoff = ENRICHMENT_PADJ,
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE
    )
  }, silent = FALSE)
  
  # Save results
  prefix <- paste0(contrast_name, "_", direction)
  
  if (!is.null(ora_go) && nrow(as.data.frame(ora_go)) > 0) {
    write.csv(as.data.frame(ora_go), 
             file.path(output_path, paste0(prefix, "_ORA_GO.csv")), 
             row.names = FALSE)
    p <- dotplot(ora_go, showCategory = 20, font.size = 8, 
                title = paste("GO BP:", contrast_name, direction))
    ggsave(file.path(output_path, paste0(prefix, "_ORA_GO_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(ora_reactome) && nrow(as.data.frame(ora_reactome)) > 0) {
    write.csv(as.data.frame(ora_reactome), 
             file.path(output_path, paste0(prefix, "_ORA_Reactome.csv")), 
             row.names = FALSE)
    p <- dotplot(ora_reactome, showCategory = 20, font.size = 8,
                title = paste("Reactome:", contrast_name, direction))
    ggsave(file.path(output_path, paste0(prefix, "_ORA_Reactome_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(ora_kegg) && nrow(as.data.frame(ora_kegg)) > 0) {
    write.csv(as.data.frame(ora_kegg), 
             file.path(output_path, paste0(prefix, "_ORA_KEGG.csv")), 
             row.names = FALSE)
    p <- dotplot(ora_kegg, showCategory = 20, font.size = 8,
                title = paste("KEGG:", contrast_name, direction))
    ggsave(file.path(output_path, paste0(prefix, "_ORA_KEGG_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  return(list(GO_BP = ora_go, Reactome = ora_reactome, KEGG = ora_kegg))
}

# Run ORA for each contrast
ora_results <- list()

for (contrast_name in names(res_list)) {
  cat(sprintf("\n=== ORA for: %s ===\n", contrast_name))
  
  sig_genes <- res_list[[contrast_name]]$sig
  
  if (nrow(sig_genes) == 0) {
    message(sprintf("No significant genes for %s", contrast_name))
    next
  }
  
  # Split by direction
  up_genes <- sig_genes %>% filter(log2FoldChange > 0) %>% pull(symbol)
  down_genes <- sig_genes %>% filter(log2FoldChange < 0) %>% pull(symbol)
  
  # Run ORA
  ora_results[[contrast_name]][["up"]] <- 
    run_ora(up_genes, contrast_name, "up", ORA_DIR)
  
  ora_results[[contrast_name]][["down"]] <- 
    run_ora(down_genes, contrast_name, "down", ORA_DIR)
}

## GENE SET ENRICHMENT ANALYSIS (GSEA) ----

cat("\n=== Running Gene Set Enrichment Analysis ===\n")

run_gsea <- function(gene_list, contrast_name, output_path) {
  
  if (length(gene_list) == 0) {
    message(sprintf("No genes for GSEA: %s", contrast_name))
    return(NULL)
  }
  
  cat(sprintf("\nRunning GSEA for %s (%d genes)\n", 
              contrast_name, length(gene_list)))
  
  # GO Biological Process
  gsea_go <- NULL
  try({
    gsea_go <- gseGO(
      geneList = gene_list,
      OrgDb = orgdb,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE,
      pvalueCutoff = ENRICHMENT_PADJ,
      eps = 0,
      seed = 123
    )
  }, silent = FALSE)
  
  # Reactome
  gsea_reactome <- NULL
  try({
    gsea_reactome <- gsePathway(
      geneList = gene_list,
      organism = ORGANISM,
      pAdjustMethod = "BH",
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE,
      pvalueCutoff = ENRICHMENT_PADJ,
      eps = 0,
      seed = 123,
      verbose = FALSE
    )
  }, silent = FALSE)
  
  # KEGG
  gsea_kegg <- NULL
  try({
    gsea_kegg <- gseKEGG(
      geneList = gene_list,
      organism = KEGG_ORGANISM,
      pAdjustMethod = "BH",
      minGSSize = MIN_GENESET_SIZE,
      maxGSSize = MAX_GENESET_SIZE,
      pvalueCutoff = ENRICHMENT_PADJ,
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
    p <- dotplot(gsea_go, showCategory = 20, font.size = 8,
                title = paste("GSEA GO:", contrast_name))
    ggsave(file.path(output_path, paste0(contrast_name, "_GSEA_GO_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
    write.csv(as.data.frame(gsea_reactome), 
             file.path(output_path, paste0(contrast_name, "_GSEA_Reactome.csv")), 
             row.names = FALSE)
    p <- dotplot(gsea_reactome, showCategory = 20, font.size = 8,
                title = paste("GSEA Reactome:", contrast_name))
    ggsave(file.path(output_path, paste0(contrast_name, "_GSEA_Reactome_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    write.csv(as.data.frame(gsea_kegg), 
             file.path(output_path, paste0(contrast_name, "_GSEA_KEGG.csv")), 
             row.names = FALSE)
    p <- dotplot(gsea_kegg, showCategory = 20, font.size = 8,
                title = paste("GSEA KEGG:", contrast_name))
    ggsave(file.path(output_path, paste0(contrast_name, "_GSEA_KEGG_dotplot.pdf")), 
          plot = p, width = 10, height = 8, dpi = 300)
  }
  
  return(list(GO_BP = gsea_go, Reactome = gsea_reactome, KEGG = gsea_kegg))
}

# Run GSEA for each contrast
gsea_results <- list()

for (contrast_name in names(res_list)) {
  cat(sprintf("\n=== GSEA for: %s ===\n", contrast_name))
  
  # Prepare ranked gene list
  gene_data <- res_list[[contrast_name]]$all %>%
    filter(!is.na(padj)) %>%
    mutate(rank = sign(log2FoldChange) * -log10(padj))
  
  # Convert to Entrez IDs
  entrez_map <- bitr(gene_data$symbol, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = orgdb)
  
  # Merge and create ranked list
  gene_data <- gene_data %>%
    inner_join(entrez_map, by = c("symbol" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    arrange(desc(rank))
  
  gene_list <- setNames(gene_data$rank, gene_data$ENTREZID)
  
  # Run GSEA
  gsea_results[[contrast_name]] <- run_gsea(gene_list, contrast_name, GSEA_DIR)
}

## SAVE WORKSPACE ----

save(ora_results, gsea_results, 
     file = file.path(GSEA_DIR, "enrichment_results.RData"))

cat("\n=== Functional Enrichment Analysis Complete ===\n")
cat(sprintf("ORA outputs: %s\n", ORA_DIR))
cat(sprintf("GSEA outputs: %s\n", GSEA_DIR))
