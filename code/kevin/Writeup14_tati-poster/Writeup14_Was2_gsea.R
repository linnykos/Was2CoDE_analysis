rm(list=ls())

library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup14/"
df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Prater_dataset_ingredients.csv")
rownames(df) <- df$Gene
df <- df[order(df$Was2_pval, decreasing = FALSE),]

Was2_logFC <- df[,"Was2_logFC"]
names(Was2_logFC) <- rownames(df)

teststat_vec <- sort(Was2_logFC, decreasing = TRUE)

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

gse_df <- as.data.frame(gse)
grep("ANXA2", gse_df$core_enrichment)

########################################

# to find which pathways to use, run GSEA on all 3 methods on Katie's data
teststat_vec <- df[,"DESeq2_logFC"]
names(teststat_vec) <- rownames(df)
teststat_vec <- sort(teststat_vec, decreasing = TRUE)
set.seed(10)
gse_deseq2 <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

teststat_vec <- df[,"eSVD_logFC"]
names(teststat_vec) <- rownames(df)
teststat_vec <- sort(teststat_vec, decreasing = TRUE)
set.seed(10)
gse_esvd <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

teststat_vec <- df[,"NEBULA_logFC"]
names(teststat_vec) <- rownames(df)
teststat_vec <- sort(teststat_vec, decreasing = TRUE)
set.seed(10)
gse_nebula <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)


########################################

# among the pathways in gse, pick the ones not in the others
pathway_others <- unique(c(gse_deseq2$ID,
                           gse_esvd$ID,
                           gse_nebula$ID))
gse_subset <- gse[!gse$ID %in%pathway_others]
gse_subset$Description
grep("ANXA2", gse_subset$core_enrichment)

selected_pathways <- gse_subset[c(7,27,47,64),"ID"]

# Filter the gseaResult object based on selected pathways
gse_subset2 <- gse
gse_subset2@result <- gse@result[gse@result$ID %in% selected_pathways, ]

plot1 <-  enrichplot::ridgeplot(gse_subset2)
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup14_Was2_gsea.png"),
                plot1,
                height = 5,
                width = 5)

plot1 <-  enrichplot::ridgeplot(gse_subset2) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup14_Was2_gsea_zoomin.png"),
                plot1,
                height = 7,
                width = 7)


