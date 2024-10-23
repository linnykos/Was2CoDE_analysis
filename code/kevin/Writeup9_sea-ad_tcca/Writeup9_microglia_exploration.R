# https://linnykos.github.io/tiltedCCA/articles/embryo.html
rm(list=ls())

library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_multiome.RData")

seurat_obj
summary(seurat_obj@meta.data)

tab_vec <- table(seurat_obj$donor_id)
stats::quantile(tab_vec)
length(unique(seurat_obj$donor_id))

Seurat::DefaultAssay(seurat_obj) <- "RNA"
length(SeuratObject::Features(seurat_obj))
length(Seurat::VariableFeatures(seurat_obj))

Seurat::DefaultAssay(seurat_obj) <- "ATAC"
length(SeuratObject::Features(seurat_obj))
length(Seurat::VariableFeatures(seurat_obj))