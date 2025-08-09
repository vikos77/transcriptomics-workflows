#################################################
# RNA-seq Data Preprocessing Script (FPKM-aware)
# Project: RmpA Regulation Analysis
# Author: Vigneshwaran Muthuraman
# Version: 2.0 - Improved FPKM handling
#################################################

# Load required libraries
library(tidyverse)
library(limma)
library(edgeR)
library(pheatmap)
library(GEOquery)
library(readxl)
library(R.utils)
library(here)
library(viridis)
library(ggrepel)

# Create logging function
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("Starting improved RmpA RNA-seq data preprocessing")

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

###################
# Data Download and Import
###################

log_message("Downloading data from GEO (GSE286114)")

# Check if data already exists to avoid redownloading
if (!file.exists("data/raw/GSE286114_family.soft.gz")) {
  tryCatch({
    gse <- getGEO("GSE286114", GSEMatrix = TRUE, destdir = "data/raw")
    metadata <- pData(phenoData(gse[[1]]))
    
    # Get supplementary files
    getGEOSuppFiles("GSE286114", baseDir = "data/raw")
    log_message("Data download successful")
  }, error = function(e) {
    log_message(paste("ERROR: Data download failed:", e$message))
    stop("Data download failed. Please check your internet connection.")
  })
} else {
  log_message("GEO data already exists, using cached version")
  gse <- getGEO("GSE286114", GSEMatrix = TRUE, destdir = "data/raw")
  metadata <- pData(phenoData(gse[[1]]))
}

# Import FPKM data
log_message("Importing FPKM expression data")

temp_file <- tempfile(fileext = ".xls")
raw_fpkm_file_gz <- "data/raw/GSE286114/GSE286114_gene_fpkm_RmpA_Escherichia.xls.gz"

if (file.exists(raw_fpkm_file_gz)) {
  tryCatch({
    R.utils::gunzip(raw_fpkm_file_gz, destname = temp_file, remove = FALSE)
    fpkm_data <- read_excel(temp_file)
    log_message(paste("Successfully imported FPKM data with", nrow(fpkm_data), "genes"))
  }, error = function(e) {
    log_message(paste("ERROR: File import failed:", e$message))
    stop("File import failed. Please check file integrity.")
  })
} else {
  stop("File not found: ", raw_fpkm_file_gz)
}

# Clean and prepare FPKM data
fpkm_matrix <- fpkm_data %>%
  select(gene_id, PLB1K1, PLB1K2, PLB1K3, `PLB1K-rmpA1`, `PLB1K-rmpA2`, `PLB1K-rmpA3`) %>%
  mutate(across(-gene_id, as.numeric)) %>%
  column_to_rownames("gene_id")

# Create sample metadata
sample_metadata <- data.frame(
  sample_id = c("PLB1K1", "PLB1K2", "PLB1K3", "PLB1K-rmpA1", "PLB1K-rmpA2", "PLB1K-rmpA3"),
  condition = factor(c(rep("control", 3), rep("rmpA_overexpression", 3)), 
                     levels = c("control", "rmpA_overexpression")),
  replicate = factor(rep(1:3, 2)),
  batch = "batch1",
  row.names = c("PLB1K1", "PLB1K2", "PLB1K3", "PLB1K-rmpA1", "PLB1K-rmpA2", "PLB1K-rmpA3")
)

# Save metadata
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")
write.csv(sample_metadata, "results/tables/sample_metadata.csv")

# Save gene information
gene_info <- fpkm_data %>%
  select(gene_id, gene_name, gene_chr, gene_start, gene_end, 
         gene_strand, gene_length, gene_biotype, gene_description)

write_csv(gene_info, "data/processed/gene_info.csv")
log_message("Saved gene annotation information")

###################
# FPKM Data Assessment and Filtering
###################

log_message("Assessing FPKM data distribution")

# Convert FPKM to TPM (better for between-sample comparisons)
fpkm_to_tpm <- function(fpkm) {
  # TPM = (FPKM / sum(FPKM)) * 10^6
  apply(fpkm, 2, function(x) (x / sum(x)) * 1e6)
}

tpm_matrix <- fpkm_to_tpm(fpkm_matrix)

# Log transformation for visualization (add pseudocount)
log_fpkm <- log2(fpkm_matrix + 0.5)
log_tpm <- log2(tpm_matrix + 0.5)

# Create distribution plots
pdf("results/figures/qc/fpkm_tpm_distributions.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

# FPKM distribution
hist(as.matrix(log_fpkm), 
     breaks = 50,
     main = "Distribution of log2(FPKM + 0.5)",
     xlab = "log2(FPKM + 0.5)",
     col = "lightblue",
     border = "white")

# TPM distribution
hist(as.matrix(log_tpm), 
     breaks = 50,
     main = "Distribution of log2(TPM + 0.5)",
     xlab = "log2(TPM + 0.5)",
     col = "lightgreen",
     border = "white")

# Sample-wise boxplots
boxplot(log_fpkm,
        main = "log2(FPKM + 0.5) per Sample",
        ylab = "log2(FPKM + 0.5)",
        las = 2,
        col = rep(c("lightblue", "lightcoral"), each = 3),
        names = gsub("PLB1K-", "R", colnames(log_fpkm)))

boxplot(log_tpm,
        main = "log2(TPM + 0.5) per Sample",
        ylab = "log2(TPM + 0.5)",
        las = 2,
        col = rep(c("lightgreen", "lightsalmon"), each = 3),
        names = gsub("PLB1K-", "R", colnames(log_tpm)))

dev.off()
log_message("Generated FPKM/TPM distribution plots")

###################
# Filtering Low Expression Genes
###################

log_message("Filtering low expression genes")

# For FPKM/TPM data, filter based on minimum expression threshold
# Genes should have FPKM > 1 in at least 3 samples (or TPM > 1)
min_fpkm <- 1
min_samples <- 3

# Apply filtering
keep_fpkm <- rowSums(fpkm_matrix >= min_fpkm) >= min_samples
keep_tpm <- rowSums(tpm_matrix >= min_fpkm) >= min_samples

# Use the intersection for conservative filtering
keep <- keep_fpkm & keep_tpm

fpkm_filtered <- fpkm_matrix[keep, ]
tpm_filtered <- tpm_matrix[keep, ]

# Log filtering results
filtering_stats <- data.frame(
  total_genes = nrow(fpkm_matrix),
  genes_kept = sum(keep),
  genes_filtered = nrow(fpkm_matrix) - sum(keep),
  percentage_filtered = round((nrow(fpkm_matrix) - sum(keep))/nrow(fpkm_matrix) * 100, 2)
)

write.csv(filtering_stats, "results/tables/filtering_statistics.csv", row.names = FALSE)

log_message(paste("Filtering results:",
                  "Total genes:", filtering_stats$total_genes,
                  "Genes kept:", filtering_stats$genes_kept,
                  "Genes filtered out:", filtering_stats$genes_filtered,
                  paste0("(", filtering_stats$percentage_filtered, "%)")))

###################
# Prepare Data for limma-voom Analysis
###################

log_message("Preparing data for limma-voom analysis")

# For limma-voom with FPKM data, we'll use log-transformed values
# Create DGEList object with FPKM values
dge_fpkm <- DGEList(counts = fpkm_filtered,
                    samples = sample_metadata,
                    genes = gene_info[match(rownames(fpkm_filtered), gene_info$gene_id), ])

# Since we have FPKM (already normalized for library size and gene length),
# we set lib.size to 1 for all samples
dge_fpkm$samples$lib.size <- rep(1e6, ncol(fpkm_filtered))

# Calculate normalization factors (TMM) to account for compositional differences
dge_fpkm <- calcNormFactors(dge_fpkm, method = "TMM")

# Save the DGEList object
saveRDS(dge_fpkm, "data/processed/dge_fpkm_object.rds")

###################
# Alternative: Prepare TPM Data
###################

# Create DGEList object with TPM values
dge_tpm <- DGEList(counts = tpm_filtered,
                   samples = sample_metadata,
                   genes = gene_info[match(rownames(tpm_filtered), gene_info$gene_id), ])

dge_tpm$samples$lib.size <- rep(1e6, ncol(tpm_filtered))
dge_tpm <- calcNormFactors(dge_tpm, method = "TMM")

saveRDS(dge_tpm, "data/processed/dge_tpm_object.rds")

###################
# Quality Control Plots
###################

log_message("Generating quality control plots")

# MDS plot
pdf("results/figures/qc/mds_plots.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# MDS for FPKM
plotMDS(dge_fpkm, 
        col = rep(c("blue", "red"), each = 3),
        main = "MDS Plot - FPKM Data",
        labels = gsub("PLB1K-", "R", colnames(fpkm_filtered)))
legend("topright", 
       legend = c("Control", "rmpA OE"), 
       col = c("blue", "red"), 
       pch = 16)

# MDS for TPM
plotMDS(dge_tpm, 
        col = rep(c("blue", "red"), each = 3),
        main = "MDS Plot - TPM Data",
        labels = gsub("PLB1K-", "R", colnames(tpm_filtered)))
legend("topright", 
       legend = c("Control", "rmpA OE"), 
       col = c("blue", "red"), 
       pch = 16)

dev.off()

# Density plots
pdf("results/figures/qc/density_plots.pdf", width = 10, height = 6)

# Create density plot function
plot_density <- function(data, title) {
  log_data <- log2(data + 0.5)
  plot(density(log_data[,1]), 
       main = title,
       xlab = "log2(Expression + 0.5)",
       ylab = "Density",
       col = 1,
       lwd = 2,
       xlim = range(log_data))
  
  for(i in 2:ncol(log_data)) {
    lines(density(log_data[,i]), col = i, lwd = 2)
  }
  
  legend("topright", 
         legend = colnames(data),
         col = 1:ncol(data),
         lty = 1,
         lwd = 2,
         cex = 0.8)
}

plot_density(fpkm_filtered, "Density Plot - Filtered FPKM Data")
dev.off()

###################
# Save Processed Data Summary
###################

# Create summary report
summary_stats <- data.frame(
  Metric = c("Total genes in dataset",
             "Genes after filtering",
             "Genes removed",
             "Percentage removed",
             "Min FPKM threshold",
             "Min samples required",
             "Median library size (FPKM sum)",
             "Data type available"),
  Value = c(nrow(fpkm_matrix),
            nrow(fpkm_filtered),
            nrow(fpkm_matrix) - nrow(fpkm_filtered),
            paste0(round((nrow(fpkm_matrix) - nrow(fpkm_filtered))/nrow(fpkm_matrix) * 100, 2), "%"),
            min_fpkm,
            min_samples,
            round(median(colSums(fpkm_matrix))),
            "FPKM (normalized)")
)

write.csv(summary_stats, "results/tables/preprocessing_summary.csv", row.names = FALSE)

# Save session info
writeLines(capture.output(sessionInfo()), "results/session_info_preprocessing.txt")

###################
# Warning and Recommendations
###################

log_message("=== IMPORTANT NOTICE ===")
log_message("This dataset contains FPKM values, which are already normalized.")
log_message("FPKM data is NOT ideal for differential expression analysis.")
log_message("")
log_message("Recommendations:")
log_message("1. Use limma-voom for differential expression (next script)")
log_message("2. Consider obtaining raw count data if possible")
log_message("3. TPM values have been calculated as an alternative")
log_message("4. Results should be interpreted with caution")
log_message("")
log_message("For robust DE analysis, raw counts are strongly preferred!")
log_message("=======================")

log_message("Preprocessing complete! Ready for limma-voom analysis.")