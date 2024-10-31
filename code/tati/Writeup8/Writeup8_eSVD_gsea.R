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


load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_esvd.RData")
esvd_res_prater <- eSVD_obj
load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_eSVD.RData")
esvd_res_seaad <- eSVD_obj
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_eSVD.RData")
esvd_res_rosmap <- eSVD_obj

intersect_genes <- Reduce(intersect, list(
  names(esvd_res_prater$teststat_vec),
  names(esvd_res_seaad$teststat_vec),
  names(esvd_res_rosmap$teststat_vec)
))

prater_df <- data.frame(
  gene = names(esvd_res_prater$case_mean),
  logFC_prater = log2(esvd_res_prater$case_mean / esvd_res_prater$control_mean)
) %>%
  filter(gene %in% intersect_genes)

seaad_df <- data.frame(
  gene = names(esvd_res_seaad$case_mean),
  logFC_seaad = log2(esvd_res_seaad$case_mean / esvd_res_seaad$control_mean)
) %>%
  filter(gene %in% intersect_genes)

rosmap_df <- data.frame(
  gene = names(esvd_res_rosmap$case_mean),
  logFC_rosmap = log2(esvd_res_rosmap$case_mean / esvd_res_rosmap$control_mean)
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
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

head(as.data.frame(gse))
gse_df <- as.data.frame(gse)
gse_df <- gse_df[gse_df$p.adjust <= 0.05,]
gse_df[1:100,"Description"]

library(tidyverse)
write_csv(gse_df, "~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_esvd_gsea_results.csv")

#ridge_plots_averaged

library(DOSE)
library(enrichplot)
library(ggplot2)

ridge_plot <- ridgeplot(gse) +
  labs(x = "Enrichment Score",
       y = "GO Biological Process",
       title = "Gene Set Enrichment Analysis Average logFC Across Datasets",
       subtitle = "Distribution of enrichment scores for top GO terms") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 8))

save_path <- "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8/Writeup8_esvd_gsea_ridge_plot.png"
ggsave(save_path, ridge_plot, width = 12, height = 8, dpi = 300)

#dot_plot
dot_plot <- enrichplot::dotplot(gse, showCategory=30) + 
  ggplot2::ggtitle("Dotplot for GSEA (Average logFC Across Datasets)")
dot_plot <- dot_plot + 
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

ggplot2::ggsave(
  filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8/Writeup8_esvd_gsea_dot_plot.png",
  dot_plot, 
  device = "png", 
  width = 7, 
  height = 7, 
  units = "in"
)


