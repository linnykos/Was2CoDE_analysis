rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(arrow)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

# Reading the Feather file
df <- arrow::read_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scvi.feather")
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

#######################

var_vec <- c("seurat_clusters", "Pt_ID", "CognitiveStatus", 
             "Study_Designation", "Sex",
             "PMI", "BrainPh", "Race", "FreshBrainWeight",
             "NIA_AA", "ThalPhase", 
             "BraakStage", "CERAD", "LATEScore", "coded_Age",
             "SeqBatch")

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup10/Writeup10_prater_scvi_umap-covariates.pdf",
    onefile = T, width = 5, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(ss_data_norm, 
                           reduction = "umap",
                           group.by = variable,
                           raster = TRUE)
  print(plot1)
}

dev.off()

##############################

# fixing all the defficiencies
# update the RNA assay
ss_data_norm[["RNA"]] <- as(object = ss_data_norm[["RNA"]], Class = "Assay5")

# change 1 in "Braak" to "BraakI", and do the same for others
# CERAD is "Absent", "Sparse", "Frequent", "Moderate"


# put factors



# add in all the highly variable genes and other relevant genes

