rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(arrow)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

# Reading the Feather file
df <- arrow::read_feather('/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-denoised.feather')
rowname_vec <- as.matrix(df[,ncol(df)])[,1]
df <- df[,-ncol(df)]
df <- as.matrix(df)
rownames(df) <- rowname_vec

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL

# remove "scale.data" slot
SeuratObject::LayerData(object = ss_data_norm, 
                        assay = "RNA", 
                        layer = "scale.data") <- NULL

Seurat::VariableFeatures(ss_data_norm) <- colnames(df)
ss_data_norm <- subset(ss_data_norm, features = Seurat::VariableFeatures(ss_data_norm))
SeuratObject::LayerData(object = ss_data_norm, 
                        assay = "RNA", 
                        layer = "data") <- t(df)

# compute PCA and UMAP
ss_data_norm <- Seurat::ScaleData(ss_data_norm)
ss_data_norm <- Seurat::RunPCA(ss_data_norm, 
                               features = Seurat::VariableFeatures(ss_data_norm),
                               verbose = FALSE)
set.seed(10)
ss_data_norm <- Seurat::RunUMAP(ss_data_norm, dims = 1:30)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- "This is Katie's data but putting the scVI-normalized data (only 2000 genes) inside under 'data'."

save(session_info, date_of_run, note,
     ss_data_norm, 
     file = "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData")


