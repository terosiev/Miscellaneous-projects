################################################################################
# Master Script: Run Complete RNA-seq Analysis Pipeline
#
# This script executes all analysis steps in sequence:
# 1. Quality control and filtering
# 2. Differential expression analysis
# 3. Functional enrichment analysis
#
# Usage: Rscript run_pipeline.R
################################################################################

# Set working directory to pipeline root
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cat("\n")
cat("================================================================================\n")
cat("               RNA-seq Differential Expression Analysis Pipeline               \n")
cat("================================================================================\n")
cat("\n")

# Track start time
start_time <- Sys.time()

# Run pipeline steps
tryCatch({
  
  cat("\n### STEP 1: Quality Control and Filtering ###\n\n")
  source("scripts/01_qc_filtering.R")
  
  cat("\n### STEP 2: Differential Expression Analysis ###\n\n")
  source("scripts/02_differential_expression.R")
  
  cat("\n### STEP 3: Functional Enrichment Analysis ###\n\n")
  source("scripts/03_functional_enrichment.R")
  
  # Calculate runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  cat("\n")
  cat("================================================================================\n")
  cat("                           PIPELINE COMPLETE                                   \n")
  cat("================================================================================\n")
  cat(sprintf("\nTotal runtime: %.2f minutes\n", as.numeric(runtime)))
  cat(sprintf("Results saved to: %s\n\n", "results/"))
  
}, error = function(e) {
  cat("\n")
  cat("================================================================================\n")
  cat("                              ERROR OCCURRED                                   \n")
  cat("================================================================================\n")
  cat("\nError message:\n")
  cat(conditionMessage(e), "\n\n")
  stop("Pipeline failed. Please check the error message above.")
})
