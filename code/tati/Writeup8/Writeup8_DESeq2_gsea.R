rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(IdeasCustom)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_Pseudobulk-DEseq2.RData")
deseq2_res_prater <- deseq2_res
load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_deseq2.RData")
deseq2_res_seaad <- deseq2_res
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_Pseudobulk-DEseq2.RData")
deseq2_res_rosmap <- deseq2_res

intersect_genes <- Reduce(intersect, list(
  rownames(deseq2_res_prater),
  rownames(deseq2_res_seaad),
  rownames(deseq2_res_rosmap)
))

prater_df <- data.frame(
  gene = rownames(deseq2_res_prater),
  logFC_prater = deseq2_res_prater[,"log2FoldChange"]
) %>%
  filter(gene %in% intersect_genes)

seaad_df <- data.frame(
  gene = rownames(deseq2_res_seaad),
  logFC_seaad = deseq2_res_seaad[,"log2FoldChange"]
) %>%
  filter(gene %in% intersect_genes)

rosmap_df <- data.frame(
  gene = rownames(deseq2_res_rosmap),
  logFC_rosmap = deseq2_res_rosmap[,"log2FoldChange"]
) %>%
  filter(gene %in% intersect_genes)

combined_logfc <- prater_df %>%
  left_join(seaad_df, by = "gene") %>%
  left_join(rosmap_df, by = "gene") %>%
  mutate(
    avg_logFC = (logFC_prater + logFC_seaad + logFC_rosmap) / 3
  ) %>%
  arrange(desc(abs(avg_logFC)))  # Sort by absolute average logFC

teststat_vec <- combined_logfc$avg_logFC
names(teststat_vec) <- combined_logfc$gene
teststat_vec <- sort(teststat_vec, decreasing = TRUE) # Sort in decreasing order

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 1,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)
head(as.data.frame(gse))
gse_df <- as.data.frame(gse)
# gse_df <- gse_df[gse_df$p.adjust <= 0.05,]
# gse_df[1:100,"Description"]

deseq2_df <- gse_df
library(tidyverse)
write_csv(deseq2_df, "~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_deseq2_gsea_results.csv")
