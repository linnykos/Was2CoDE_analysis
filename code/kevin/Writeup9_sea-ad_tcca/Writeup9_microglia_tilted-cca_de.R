rm(list=ls())

library(Seurat)
library(tiltedCCA)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_tcca-with-synchrony.RData")

# https://satijalab.org/seurat/articles/de_vignette
# Find DE features between CD16 Mono and CD1 Mono
Idents(seurat_obj) <- "RNA_snn_res.0.5"
sapply(levels(Idents(seurat_obj)), function(lev){
  mean(seurat_obj$synchrony[which(seurat_obj$RNA_snn_res.0.5==lev)])
})

set.seed(10)
Seurat::DefaultAssay(seurat_obj) <- "RNA"
synchrony_markers <- Seurat::FindMarkers(seurat_obj, 
                                         ident.1 = "11",
                                         test.use = "wilcox")

sig_genes <- rownames(synchrony_markers)[synchrony_markers$p_val_adj <= 0.05]
sig_genes_positive <- intersect(
  sig_genes,
  rownames(synchrony_markers)[synchrony_markers$avg_log2FC > 0]
)
sig_genes_negative <- intersect(
  sig_genes,
  rownames(synchrony_markers)[synchrony_markers$avg_log2FC < 0]
)
# convert gene names
edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

sig_genes2 <- AnnotationDbi::mapIds(edb,
                                    keys = sig_genes,
                                    keytype = "GENEID", 
                                    column = "GENENAME", 
                                    multiVals = "first")
names(sig_genes2) <- NULL

write.table(
  sig_genes2,
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE,
  file = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/csv/kevin/Writeup9/Writeup9_microglia_synchrony-markers.csv"
)

sig_genes_positive2 <- AnnotationDbi::mapIds(edb,
                                    keys = sig_genes_positive,
                                    keytype = "GENEID", 
                                    column = "GENENAME", 
                                    multiVals = "first")
write.table(
  sig_genes_positive2,
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE,
  file = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/csv/kevin/Writeup9/Writeup9_microglia_synchrony-markers_positive.csv"
)

sig_genes_negative2 <- AnnotationDbi::mapIds(edb,
                                             keys = sig_genes_negative,
                                             keytype = "GENEID", 
                                             column = "GENENAME", 
                                             multiVals = "first")
write.table(
  sig_genes_negative2,
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE,
  file = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/csv/kevin/Writeup9/Writeup9_microglia_synchrony-markers_negative.csv"
)

