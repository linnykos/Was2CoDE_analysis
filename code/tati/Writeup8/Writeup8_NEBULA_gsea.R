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


load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_NEBULA.RData")
nebula_res_prater <- nebula_res
load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_NEBULA.RData")
nebula_res_seaad <- nebula_res
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_NEBULA.RData")
nebula_res_rosmap <- nebula_res

intersect_genes <- Reduce(intersect, list(
  nebula_res_prater$summary$gene,
  nebula_res_seaad$summary$gene,
  nebula_res_rosmap$summary$gene
))

prater_df <- nebula_res_prater$summary %>%
  filter(gene %in% intersect_genes) %>%
  dplyr::select(gene, logFC_Study_DesignationAD) %>%
  rename(logFC_prater = logFC_Study_DesignationAD)

seaad_df <- nebula_res_seaad$summary %>%
  filter(gene %in% intersect_genes) %>%
  dplyr::select(gene, logFC_ADNCCase) %>%
  rename(logFC_seaad = logFC_ADNCCase)

rosmap_df <- nebula_res_rosmap$summary %>%
  filter(gene %in% intersect_genes) %>%
  dplyr::select(gene, logFC_ADpathyes) %>%
  rename(logFC_rosmap = logFC_ADpathyes)


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

nebula_df <- gse_df
library(tidyverse)
write_csv(nebula_df, "~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_nebula_gsea_results.csv")
