# Quality control for longitudinal RNA-seq analysis
# Author: Tero Sievänen
# Description: QC workflow including PCA, hierarchical clustering, and data filtering

# Load libraries ----
library(DESeq2)
library(edgeR)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(tximport)
library(limma)

# Define paths ----
# NOTE: Update these paths to match your directory structure
output_path_qc <- "results/qc/"
pdf_filename <- "PCA"

# Define functions ----

# Principal component analysis
run_PCA <- function(dds, design_matrix, pdf_filename = NULL, path = NULL) {
  
  # Take rlog to further scale data (required a generated dds-object)
  rld <- rlog(dds, blind = TRUE)
  
  # Perform PCA using prcomp
  pca <- prcomp(t(assay(rld)))
  
  # Combine design matrix with PCA results
  df <- cbind(design_matrix, pca$x)
  
  # Define intgroup
  intgroup <- c("batch","time", "condition")
  
  # Create separate plots for each group
  plots <- lapply(intgroup, function(group) {
    plot <- ggplot(df, aes(x = PC1, y = PC2, color = .data[[group]])) + 
      geom_point(size = 5) +
      geom_text(aes(label = rownames(df)), nudge_x = 0.4, nudge_y = 0.4, check_overlap = TRUE, color = "black") +
      theme_classic() +
      ggtitle(paste("PCA Plot (Group:", group, ")"))
    
    # Save plot as PDF if filename provided
    if (!is.null(pdf_filename)) {
      filename <- if (!is.null(path)) file.path(path, paste0(pdf_filename, "_", group, ".pdf")) else paste0(pdf_filename, "_", group, ".pdf")
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

# Data setup ----

# Import count data
# NOTE: Update path to your count matrix
counts <- read.csv("data/counts_matrix.csv", sep = "\t", row.names = 1)

# Import metadata
# NOTE: Update path to your metadata file
meta <- read.csv("data/metadata.csv", sep = "\t", row.names = 1)

# Load the annotation table for GrCh38
# NOTE: Update path to your tx2gene file
tx2gene <- read.delim("data/tx2gene.tsv", header = F)

# Create an annotation file
grch38annot <- tx2gene %>%
  dplyr::select(gene = V2, symbol = V3) %>%
  distinct() # Remove duplicate genes

# Optional: Filter out specific samples if needed
# Example: exclude outliers or high-risk samples
# select <- which(meta$sample_id == "SAMPLE_X" | meta$risk_category == "HIGH_RISK")
# meta <- meta[-select, ]
# counts <- counts[, -select]

# Extract interesting data from meta file and create an annotation file
targets <- data.frame(
  condition = as.factor(meta$condition),
  time = as.factor(meta$timepoint),
  replicate = as.factor(meta$subject_id),
  batch = as.factor(meta$batch),
  row.names = colnames(counts))

# Create a data matrix
y <- as.matrix(counts)
rownames(y) <- rownames(counts)

# Create a DGEList
data <- DGEList(y)

# Add sample information
data$samples$group <- factor(rownames(targets))
data$samples$batch <- factor(targets$batch)
data$samples$subject <- factor(targets$replicate)
levels(targets$time) <- c("time1",  "time2","time3")
data$samples$time <- factor(targets$time)
data$samples$condition <- factor(targets$condition)

# Add gene name information
data$genes <- data.frame(grch38annot$symbol[match(rownames(data), grch38annot$gene)])

# Check final data format
data



# Quality control ----

# Create DESeq2DataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = targets, design = ~ batch + time + condition)

# Compute size factors
dds <- estimateSizeFactors(dds)

# Perform quality controls and produce QC-plots
run_PCA(dds, targets, pdf_filename, output_path_qc)

# Hierarchical clustering and sample-wise correlations

# Compute rlog transformed counts
rld <- rlog(dds, blind=T)

# Extract
rld_mat <- assay(rld)

# Compute pairwise correlation values based on rlog values
rld_cor <- cor(assay(rld))

# Choose annotations
anno <- targets[, -c(3)]

# Set a color palette
heat_colors <- brewer.pal(5, "RdBu")

# Specify colors
ann_colors = list(
  time = c(time1 ="#66C2A5", 
           time2="#FC8D62",
           time3="#8DA0CB"),
  batch = c(B1 = "#E78AC3", B2 = "#A6D854"),
  condition = c(CONDITION1 = "#D53E4F", CONDITION2 = "#3288BD")
)

# Draw heatmap of sample-wise correlations
heatmap_plot <- pheatmap(rld_cor, annotation = anno, annotation_colors = ann_colors, cluster_rows = F, cluster_cols = F, angle_col = 45, fontsize = 7)
ggsave(filename = file.path(output_path_qc, "qc_hmap.pdf"), heatmap_plot, dpi = 300)


# Data filtering ----

# Transform to cpm
lcpm <- cpm(data, log=TRUE)

# The log-CPM values will be used for exploratory plots. When log=TRUE, the cpm function adds an offset to the CPM values before converting to the log2-scale. 
# By default, the offset is 2/L where 2 is the "prior count" and L is the average library size in millions, so the log-CPM values are related to the CPM values by log2(CPM + 2/L).
# This calculation ensures that any two read counts with identical CPM values will also have identical log-CPM values. 
# The prior count avoids taking the logarithm of zero, and also reduces spurious variability for genes with very low counts by shrinking all the inter-sample log-fold-changes towards zero, something that is helpful for exploratory plotting.

L <- mean(data$samples$lib.size) * 1e-6
M <- median(data$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

# Remove low expressed genes
table(rowSums(data$counts==10)==2)

paste0("There are ", sum(table(rowSums(data$counts==0)==2)[1]), " expressed genes and ", 
       sum(table(rowSums(data$counts==0)==2)[2]), " unexpressed genes in the dataset.")

paste0("Unexpressed genes account for ", 
       round(sum(table(rowSums(data$counts==0)==2)[2])/(sum(table(rowSums(data$counts==0)==2)[1])+sum(table(rowSums(data$counts==0)==2)[2]))*100, digits=2),
       "% of the genes in the data.")

# filterByExpr provides a way to filter genes, while keeping as many genes as possible with worthwhile counts.
keep.exprs <- edgeR::filterByExpr(data, min.count = 30, min.total.count = 50)
data <- data[keep.exprs,, keep.lib.sizes = F]
dim(data)


# Visualizing the filtering done by filterByExpr

# Set options
options(repr.plot.width=16, repr.plot.height=8)

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(lcpm)
col <- brewer.pal(min(nsamples, 12), "Set3")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
lcpm <- cpm(data, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

# Inspect samples with MDS plots

# Set options for plotting
par(mfrow=c(1,2))

# Convert factor levels to character strings
col.group <- as.character(data$samples$group)
col.batch <- as.character(data$samples$batch)

# Specify colors based on factor levels
set1_palette <- brewer.pal(min(nlevels(data$samples$group), 9), "Set1")
set2_palette <- brewer.pal(nlevels(data$samples$batch), "Set2")

# Map factor levels to colors
col.group <- set1_palette[as.numeric(factor(col.group))]
col.batch <- set2_palette[as.numeric(factor(col.batch))]

# Plot MDS with specified colors
plotMDS(lcpm, labels=data$samples$group, col=col.group)
title(main="A. Samples")
plotMDS(lcpm, labels=data$samples$batch, col=col.batch)
title(main="B. Batch")


# Data normalization ----

# During the sample preparation or sequencing process, external factors that are not of biological interest can affect the expression of individual samples. 
# For example, samples processed in the first batch of an experiment can have higher expression overall when compared to samples processed in a second batch. 
# It is assumed that all samples should have a similar range and distribution of expression values. Normalization is required to ensure that the expression distributions of each sample are similar across the entire experiment.
par(mfrow=c(1,2), mar=c(12, 4, 4, 2))
lcpm <- cpm(data, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis = 0.5)
title(main="A. Unnormalised data",ylab="Log-cpm")
data <- calcNormFactors(data, method = "TMM")  
data$samples$norm.factors
lcpm <- cpm(data, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis = 0.5)
title(main="B. Normalised data",ylab="Log-cpm")
# Reset margins to default
par(mar=c(5, 4, 4, 2))
