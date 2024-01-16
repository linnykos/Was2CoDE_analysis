rm(list=ls())
library(Seurat)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData")

var_vec <- c("seurat_clusters", "Pt_ID", "CognitiveStatus", 
             "Study_Designation", "Sex", "genotype_APOE",
             "PMI", "BrainPh", "Race", "FreshBrainWeight",
             "NIA_AA", "ThalPhase", 
             "BraakStage", "CERAD", "LATEScore", "coded_Age",
             "SeqBatch")

pdf("../../../figures/kevin/Writeup1/naive-preprocess_umap-covariates.pdf",
    onefile = T, width = 5, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(ss_data_norm, 
                           reduction = "umap",
                           group.by = variable,
                           raster = TRUE)
  print(plot1)
}

dev.off()