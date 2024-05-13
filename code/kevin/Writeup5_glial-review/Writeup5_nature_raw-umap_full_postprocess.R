rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_raw-umap_full.RData")

ss_data_norm[["umap"]] <- seurat_umap

var_vec <- c("seurat_clusters", "Pt_ID", "CognitiveStatus", 
             "Study_Designation", "Sex", "genotype_APOE",
             "PMI", "BrainPh", "Race", "FreshBrainWeight",
             "NIA_AA", "ThalPhase", 
             "BraakStage", "CERAD", "LATEScore", "coded_Age",
             "SeqBatch")
var_vec <- var_vec[var_vec %in% colnames(ss_data_norm@meta.data)]

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_nature_raw-umap_covariates.pdf",
    onefile = T, width = 8, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(ss_data_norm, 
                           label = TRUE,
                           repel = TRUE,
                           label.size = 2.5,
                           reduction = "umap",
                           group.by = variable,
                           raster = TRUE)
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  print(plot1)
}

dev.off()