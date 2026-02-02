################################################################################
# Package Installation Script
#
# This script installs all required R packages for the RNA-seq pipeline
################################################################################

cat("\n=== Installing Required R Packages ===\n\n")

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  cat("✓ Installed BiocManager\n")
}

# Define required packages
bioc_packages <- c(
  "DESeq2",
  "edgeR",
  "clusterProfiler",
  "ReactomePA",
  "enrichplot",
  "org.Mm.eg.db",  # Mouse (change/add others as needed)
  "apeglm",
  "ashr",
  "pheatmap"
)

cran_packages <- c(
  "tidyverse",
  "data.table",
  "ggrepel",
  "RColorBrewer",
  "gridExtra",
  "msigdbr"
)

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

# Install CRAN packages
cat("\nInstalling CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE)
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

cat("\n=== Package Installation Complete ===\n")
cat("\nFor other organisms, install the appropriate annotation package:\n")
cat("  BiocManager::install('org.Hs.eg.db')  # Human\n")
cat("  BiocManager::install('org.Dm.eg.db')  # Fly\n")
cat("  BiocManager::install('org.Rn.eg.db')  # Rat\n\n")
