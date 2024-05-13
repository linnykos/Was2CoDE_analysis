rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL

seurat_obj <- ss_data_norm
rm(ls = "ss_data_norm"); gc(TRUE)

seurat_obj$Pt_ID <- paste0("D:", as.character(seurat_obj$Pt_ID))
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[which(seurat_obj$Pt_ID %in% c("D:1", "D:2", "D:3", "D:9", 
                                       "D:21", "D:22", "D:4", "D:5", 
                                       "D:6", "D:11"))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# https://github.com/satijalab/seurat/issues/8183
# https://github.com/satijalab/seurat/issues/7501#issuecomment-1854571904
print("Performing SCTransform")
set.seed(10)
unwanted <- c("percent.mito")
nfeatures <- 2000
seurat_obj <- Seurat::SCTransform(seurat_obj, 
                                  variable.features.n = nfeatures,
                                  do.scale = FALSE,
                                  do.center = FALSE,
                                  conserve.memory = TRUE,
                                  return.only.var.genes = FALSE,
                                  vars.to.regress = unwanted,
                                  verbose = FALSE)

print("Performing PCA")
set.seed(10)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             features = Seurat::VariableFeatures(seurat_obj), 
                             npcs = 50,
                             verbose = FALSE)

print("Performing UMAP")
set.seed(10)
seurat_obj <- Seurat::RunUMAP(object = seurat_obj, dims = 1:30)
seurat_umap <- seurat_obj[["umap"]]
metadata_df <- seurat_obj@meta.data

save(seurat_umap, metadata_df,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_raw-umap.RData")




