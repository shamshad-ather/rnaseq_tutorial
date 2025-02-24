# Install BiocManager package/library
install.packages("BiocManager")

# Install Bioconductor packages/library
BiocManager::install(c("DESeq2", "EnhancedVolcano", "ComplexHeatmap"))

# Install general packages
install.packages(c("tidyr", "dplyr"))

# Load the Libraries
library(DESeq2)           # For differential expression analysis
library(tidyr)            # For data wrangling
library(dplyr)            # For data manipulation
library(EnhancedVolcano)  # For visualization of DEGs
library(ComplexHeatmap)   # For heatmap visualization

#set working directory
setwd("D:/RNAseq_tutorial/rawdata/")  # Replace with your actual path

# Check if the working directory is set correctly
getwd()

# Load the raw counts, skipping the first comment lines starting with '#'
counts_raw <- read.delim("counts.txt", row.names = 1, comment.char = "#")

# Check the first few rows
head(counts_raw)

# Select only the Geneid and BAM count columns
counts <- counts_raw[, c(6:ncol(counts_raw))]


# Check the first few rows
head(counts)

# Check sample names
colnames(counts)

# clean the column name
colnames(counts) <- gsub("_sorted.bam", "", colnames(counts))

# Rename SRR* to test and control
colnames(counts) <- c("CTRL_1", "CTRL_2", "TEST_1", "TEST_2")

# Check the updated column names
print(colnames(counts))

# Create a simple condition vector
condition <- factor(c("TEST", "TEST", "CTRL", "CTRL"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = data.frame(condition), 
                              design = ~ condition)
# dds view
dds
colData(dds)  # Check sample conditions
assay(dds)[1:5, ]  # View first 5 genes' counts

# Perform normalization and statistical analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "TEST", "CTRL"))  # Extract results

# filter DE genes
res_filtered <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# summary and count of DE genes
summary(res)
nrow(res_filtered)

# Save all results
write.csv(res_filtered, "DEGs_filtered.csv")
write.csv(as.data.frame(res), "DEGs_full.csv")

# Visualize the result
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot',
                subtitle = 'TEST vs CTRL',
                caption = 'DESeq2 results')

# Heatmap the DE genes
# Select genes with padj > 0.05
significant_genes <- res[which(res$padj < 0.05), ]

# Extract normalized counts for these genes
normalized_counts <- counts(dds, normalized=TRUE)
heatmap_data <- normalized_counts[rownames(significant_genes), ]

# Scale row-wise (z-score normalization)
heatmap_data_scaled <- t(scale(t(heatmap_data)))  # Scale across rows

Heatmap(heatmap_data_scaled, name = "Z-score", cluster_rows = TRUE, cluster_columns = TRUE,show_row_dend = F,show_column_dend = F)


