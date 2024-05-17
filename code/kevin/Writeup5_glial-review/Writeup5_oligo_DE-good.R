rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")

seurat_all <- subset(seurat_all, Supertype == "Oligo_2")

quantile(seurat_all$Continuous.Pseudo.progression.Score)
table(seurat_all$Continuous.Pseudo.progression.Score, seurat_all$donor_id)

median_val <- stats::median(seurat_all$Continuous.Pseudo.progression.Score)
donor_vec <- rep("low", length(Seurat::Cells(seurat_all)))
donor_vec[seurat_all$Continuous.Pseudo.progression.Score > median_val] <- "high"
seurat_all$classification <- as.factor(donor_vec)

set.seed(10)
Seurat::Idents(seurat_all) <- "classification"
de_markers <- Seurat::FindMarkers(seurat_all, 
                                  ident.1 = "low", 
                                  ident.2 = "high",
                                  min.pct = 0,
                                  logfc.threshold = 0)
# view results
head(de_markers)
dim(de_markers)

edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
genename_vec <- AnnotationDbi::mapIds(
  edb,
  keys = rownames(de_markers),
  keytype = "GENEID", 
  column = "GENENAME",
  multiVals = "first"
)
na_idx <- which(is.na(genename_vec))
if(length(na_idx) > 0){
  de_markers <- de_markers[-na_idx,]
  genename_vec <- genename_vec[-na_idx]
}

rownames(de_markers) <- genename_vec

# determine p-value cutoff
idx <- which(de_markers$p_val_adj <= 0.05)
pCutoff <- max(de_markers$p_val[idx])

# determine log fold change cutoff
FCcutoff <- sort(abs(de_markers$avg_log2FC), decreasing = TRUE)[200]

plot1 <- EnhancedVolcano::EnhancedVolcano(de_markers,
                                          lab = rownames(de_markers),
                                          x = "avg_log2FC",
                                          y = "p_val",
                                          pCutoff = pCutoff,
                                          FCcutoff = FCcutoff,
                                          title = 'Volcano Plot',
                                          subtitle = 'Comparison of Dementia vs. No Dementia',
                                          xlab = 'Log2 Fold Change',
                                          ylab = '-Log10 pval_ideas')

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_oligo_DE-good.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

plot1 <- EnhancedVolcano::EnhancedVolcano(de_markers,
                                          lab = rownames(de_markers),
                                          x = "avg_log2FC",
                                          y = "p_val",
                                          ylim = c(0,60),
                                          pCutoff = pCutoff,
                                          FCcutoff = FCcutoff,
                                          title = 'Volcano Plot',
                                          subtitle = 'Comparison of Dementia vs. No Dementia',
                                          xlab = 'Log2 Fold Change',
                                          ylab = '-Log10 pval_ideas')

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_oligo_DE-good_ylim.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")


