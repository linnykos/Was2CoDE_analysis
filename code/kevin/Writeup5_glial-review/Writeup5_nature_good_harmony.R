rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL

seurat_obj <- ss_data_norm
ls_vec <- ls(); ls_vec <- setdiff(ls_vec, "seurat_obj")
rm(ls = ls_vec); gc(TRUE)

seurat_obj$Pt_ID <- paste0("D:", as.character(seurat_obj$Pt_ID))
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[which(seurat_obj$Pt_ID %in% c("D:14", "D:20", "D:21", "D:18", 
                                       "D:22", "D:13", "D:4", "D:1", 
                                       "D:12", "D:11"))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

print("Splitting the data")
batch1_obj <- subset(seurat_obj, SeqBatch == 1)
batch2_obj <- subset(seurat_obj, SeqBatch == 2)

print("Splitting by donors")
batch1_obj[["RNA"]] <- split(batch1_obj[["RNA"]], f = batch1_obj$Pt_ID)
batch2_obj[["RNA"]] <- split(batch2_obj[["RNA"]], f = batch2_obj$Pt_ID)

print("Normalizing both halves")
batch1_obj <- Seurat::NormalizeData(batch1_obj)
batch1_obj <- Seurat::FindVariableFeatures(batch1_obj)
batch2_obj <- Seurat::NormalizeData(batch2_obj)
batch2_obj <- Seurat::FindVariableFeatures(batch2_obj)

print("Aggregating variable features")
gene_vec <- sort(unique(c(Seurat::VariableFeatures(batch1_obj), 
                          Seurat::VariableFeatures(batch2_obj))))
Seurat::VariableFeatures(batch1_obj) <- gene_vec
Seurat::VariableFeatures(batch2_obj) <- gene_vec

print("Running PCA on each half")
set.seed(10)
batch1_obj <- Seurat::ScaleData(batch1_obj)
batch1_obj <- Seurat::RunPCA(batch1_obj, 
                             verbose = FALSE)
set.seed(10)
batch2_obj <- Seurat::ScaleData(batch2_obj)
batch2_obj <- Seurat::RunPCA(batch2_obj, 
                             verbose = FALSE)

print("Running Harmony on each half")
batch1_obj <- Seurat::IntegrateLayers(object = batch1_obj, 
                                      method = Seurat::HarmonyIntegration, 
                                      orig.reduction = "pca",
                                      new.reduction = "harmony", 
                                      verbose = TRUE)
batch2_obj <- Seurat::IntegrateLayers(object = batch2_obj, 
                                      method = Seurat::HarmonyIntegration, 
                                      orig.reduction = "pca",
                                      new.reduction = "harmony", 
                                      verbose = TRUE)

print("Manually combine both batches into one big dataset")
harmony1 <- batch1_obj[["harmony"]]@cell.embeddings
harmony2 <- batch2_obj[["harmony"]]@cell.embeddings
harmony_all <- rbind(harmony1, harmony2)
colnames(harmony_all) <- paste0("pc_", 1:ncol(harmony_all))

metadata1 <- batch1_obj@meta.data
metadata2 <- batch2_obj@meta.data
metadata_all <- rbind(metadata1, metadata2)

seurat_all <- Seurat::CreateSeuratObject(counts = t(harmony_all), 
                                         data = t(harmony_all), 
                                         meta.data = metadata_all)
seurat_all[["RNA"]] <- split(seurat_all[["RNA"]], f = seurat_all$SeqBatch)
seurat_all[["pca"]] <- Seurat::CreateDimReducObject(embeddings = harmony_all)
Seurat::VariableFeatures(seurat_all) <- SeuratObject::Features(seurat_all)
seurat_all <- Seurat::ScaleData(seurat_all)

print("Running one last Harmony")
seurat_all <- Seurat::IntegrateLayers(object = seurat_all, 
                                      method = Seurat::HarmonyIntegration, 
                                      orig.reduction = "pca",
                                      new.reduction = "harmony", 
                                      verbose = TRUE,
                                      features = SeuratObject::Features(seurat_all))

set.seed(10)
seurat_all <- Seurat::RunUMAP(object = seurat_all, 
                              reduction = "harmony",
                              dims = 1:30)


var_vec <- c("Pt_ID", "Study_Designation", "SeqBatch")

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(seurat_all, 
                           reduction = "umap",
                           label = TRUE,
                           repel = TRUE,
                           label.size = 2.5,
                           group.by = variable)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("UMAP, Good design, after Harmony"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  ggplot2::ggsave(plot1, 
                  file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_nature_good_harmony_umap-", variable, ".png"),
                  height = 5, width = 8,
                  units = "in")
}



