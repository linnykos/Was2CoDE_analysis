load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_NEBULA.RData")
head(nebula_res$summary)

# teststat_vec <- -log10(nebula_res$summary[,"p_Study_DesignationAD"]) # our "proxy" for a test-statistics
teststat_vec <- nebula_res$summary[,"logFC_ADpathyes"]
names(teststat_vec) <- nebula_res$summary[,"gene"]

# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("ggtree")
# BiocManager::install("enrichplot")
# BiocManager::install("clusterProfiler")
# BiocManager::install("AnnotationDbi")

library(org.Hs.eg.db)
library(clusterProfiler)

teststat_vec <- sort(teststat_vec, decreasing = TRUE)

# # Convert symbols to check if they map correctly to gene IDs
# mapped_genes <- clusterProfiler::bitr(names(teststat_vec), 
#                                       fromType = "SYMBOL", 
#                                       toType = "ENTREZID", 
#                                       OrgDb = org.Hs.eg.db)
# 
# # Filter teststat_vec to keep only the mapped genes
# teststat_vec_mapped <- teststat_vec[names(teststat_vec) %in% mapped_genes$SYMBOL]
# 
# # Replace gene symbols with ENTREZ IDs
# names(teststat_vec_mapped) <- mapped_genes$ENTREZID[match(names(teststat_vec_mapped), mapped_genes$SYMBOL)]
# 
# # Sort the updated vector (optional if needed)
# teststat_vec_mapped <- sort(teststat_vec_mapped, decreasing = TRUE)
# head(teststat_vec_mapped)

# https://www.youtube.com/watch?v=Mi6u4r0lJvo
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

plot1 <- enrichplot::dotplot(gse, showCategory=30) + ggplot2::ggtitle("dotplot for GSEA (ROSMAP&NEBULA)")
plot1 <- plot1 + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/Writeup7_GSEA_ROSMAP_NEBULA_dotplot.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html?q=gsea#ridgeline-plot-for-expression-distribution-of-gsea-result
# clusterProfiler::dotplot(gse, showCategory=30) + ggtitle("dotplot for GSEA")

load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_eSVD.RData")
head(eSVD_obj)
teststat_vec <- eSVD_obj[,log2(eSVD_obj$case_mean / eSVD_obj$control)]
names(teststat_vec) <- eSVD_obj[,"gene"]
