# scripts/dge_analysis.R

# Load libraries
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Set directory path for input files
parent_dir <- "E:/GDC_TOOL/GDC_Downloads"

#--------------------------#
# STEP 1: COUNT MATRIX
#--------------------------#

tsv_files <- list.files(path = parent_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

count_list <- lapply(tsv_files, function(file) {
  df <- read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    filter(str_starts(gene_id, "ENSG")) %>%
    mutate(gene_id = str_replace(gene_id, "\\.\\d+$", "")) %>%
    select(GeneID = gene_id, Count = stranded_first)
  
  sample_name <- basename(dirname(file))
  df <- rename(df, !!sample_name := Count)
  return(df)
})

count_matrix <- reduce(count_list, full_join, by = "GeneID")
count_matrix[is.na(count_matrix)] <- 0
write_tsv(count_matrix, "cleaned_count_matrix.tsv")

#--------------------------#
# STEP 2: METADATA
#--------------------------#

metadata <- data.frame(
  SampleID = colnames(count_matrix)[-1],
  Condition = c(rep("Tumor", 10), rep("Normal", 10))
)
write.csv(metadata, "sample_metadata.csv", row.names = FALSE)

#--------------------------#
# STEP 3: DESEQ2 ANALYSIS
#--------------------------#

count_matrix <- read_tsv("cleaned_count_matrix.tsv")
metadata <- read.csv("sample_metadata.csv")

dds <- DESeqDataSetFromMatrix(countData = count_matrix[,-1], 
                              colData = metadata, 
                              design = ~ Condition)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "Tumor", "Normal"))

# Save significant DEGs
res <- res[!is.na(res$padj), ]
sig_DEGs <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
sig_DEGs$gene <- rownames(sig_DEGs)
write.csv(sig_DEGs, "significant_DEGs.csv", row.names = FALSE)

#--------------------------#
# STEP 4: PLOTS
#--------------------------#

# Volcano Plot
volcano_data <- data.frame(log2FoldChange = res$log2FoldChange,
                           pvalue = res$pvalue)
volcano_data$significant <- ifelse(volcano_data$pvalue < 0.05 & abs(volcano_data$log2FoldChange) > 1, "yes", "no")

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  scale_color_manual(values = c("gray", "red"))

# MA Plot
plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))

# Heatmap
top_genes <- head(order(res$padj), 20)
pheatmap(assay(dds)[top_genes, ], cluster_rows = TRUE, cluster_cols = TRUE)

# PCA Plot
vst_data <- vst(dds, blind = FALSE)
plotPCA(vst_data, intgroup = "Condition")

#--------------------------#
# STEP 5: FILTER TOP DEGs
#--------------------------#

top_upregulated <- res %>%
  as.data.frame() %>%
  filter(baseMean > 10) %>%
  arrange(desc(log2FoldChange)) %>%
  head(100)
top_upregulated$ensembl_gene_id <- rownames(top_upregulated)
write.csv(top_upregulated, "top_upregulated.csv", row.names = FALSE)

#--------------------------#
# STEP 6: SURFACE GENES
#--------------------------#

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
surface_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "go", values = "GO:0005886", mart = ensembl
)

surface_degs <- top_upregulated %>%
  filter(ensembl_gene_id %in% surface_genes$ensembl_gene_id) %>%
  left_join(surface_genes, by = "ensembl_gene_id")
write.csv(surface_degs, "surface_degs.csv", row.names = FALSE)

#--------------------------#
# STEP 7: CANDIDATE GENE EXPRESSION
#--------------------------#

target_genes <- c("CEACAM6", "ERBB2", "MUC1")
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$ensembl_gene_id <- rownames(norm_counts_df)

gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = rownames(norm_counts),
                  mart = ensembl)

norm_counts_named <- norm_counts_df %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%
  filter(external_gene_name %in% target_genes)

long_counts <- norm_counts_named %>%
  pivot_longer(cols = -c(ensembl_gene_id, external_gene_name),
               names_to = "sample", values_to = "expression")

sample_metadata <- as.data.frame(colData(dds))
sample_metadata$sample <- rownames(sample_metadata)
colnames(sample_metadata) <- tolower(colnames(sample_metadata))

long_counts <- long_counts %>%
  left_join(sample_metadata[, c("sample", "condition")], by = "sample")

ggplot(long_counts, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ external_gene_name, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Expression of Candidate CAR-T Target Genes",
       x = "Condition", y = "Normalized Expression") +
  scale_fill_manual(values = c("Tumor" = "#E64B35", "Normal" = "#4DBBD5"))
