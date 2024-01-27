rm(list=ls())
library("Seurat")

seurat_obj <- readRDS("~/kzlinlab/data/sea-ad/microglia-pvm_dpc.rds")
seurat_obj

# do a minor adjustment
covariate_vec <- as.character(seurat_obj$`LATE-NC stage`)
idx <- which(covariate_vec == "Staging Precluded by FTLD with TDP43 or ALS/MND or TDP-43 pathology is unclassifiable")
covariate_vec[idx] <- "Unclassifiable"
seurat_obj$`LATE-NC stage` <- factor(covariate_vec)

# do some housekeeping
seurat_obj[["scVI"]] <- NULL
seurat_obj[["umap"]] <- NULL

Seurat::DefaultAssay(seurat_obj) <- "RNA"
set.seed(10)
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                              selection.method = "vst",
                                              nfeatures = 2000)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj,
                             features = Seurat::VariableFeatures(object = seurat_obj),
                             verbose = F)
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              dims = 1:30)

var_vec <- c("Age at death", "Cognitive status",
             "Braak stage", "Thal phase",
             "CERAD score", "APOE4 status",
             "LATE-NC stage")

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup2/Writeup2_microglia-pvm_eda_umap.pdf",
    onefile = T, width = 5, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = "umap",
                           group.by = variable)
  print(plot1)
}

dev.off()

print("Done! :)")
