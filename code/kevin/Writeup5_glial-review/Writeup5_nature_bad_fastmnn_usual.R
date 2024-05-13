rm(list=ls())
library(Seurat)
library(SeuratWrappers)
library(batchelor)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL

seurat_obj <- ss_data_norm
ls_vec <- ls(); ls_vec <- setdiff(ls_vec, "seurat_obj")
rm(ls = ls_vec); gc(TRUE)

keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[which(seurat_obj$Pt_ID %in% c("D:14", "D:20", "D:21", "D:18", 
                                       "D:19", "D:15", "D:7", "D:10", 
                                       "D:4", "D:12"))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

print("Splitting by donors")
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$Pt_ID)

print("Normalizing both halves")
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)

print("Running PCA on each half")
set.seed(10)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             verbose = FALSE)

print("Running fastMNN on each half")
seurat_obj <- Seurat::IntegrateLayers(object = seurat_obj, 
                                      method = SeuratWrappers::FastMNNIntegration, 
                                      orig.reduction = "pca",
                                      new.reduction = "fastmnn", 
                                      verbose = TRUE)

set.seed(10)
seurat_obj <- Seurat::RunUMAP(object = seurat_obj, 
                              reduction = "fastmnn",
                              dims = 1:30)

var_vec <- c("Pt_ID", "Study_Designation", "SeqBatch")

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = "umap",
                           label = TRUE,
                           repel = TRUE,
                           label.size = 2.5,
                           group.by = variable)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("UMAP, Bad design, after FastMNN (usual)"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  ggplot2::ggsave(plot1, 
                  file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_nature_bad_fastmnn-usual_umap-", variable, ".png"),
                  height = 5, width = 8,
                  units = "in")
}




