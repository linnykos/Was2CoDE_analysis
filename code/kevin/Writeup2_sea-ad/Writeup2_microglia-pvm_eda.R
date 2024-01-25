rm(list=ls())
library("Seurat")

seurat_obj <- readRDS("~/kzlinlab/data/sea-ad/microglia-pvm_dpc.rds")
seurat_obj

# do a minor adjustment
covariate_vec <- seurat_obj$`LATE-NC stage`
level_vec <- levels(covariate_vec)
idx <- which(level_vec == "Staging Precluded by FTLD with TDP43 or ALS/MND or TDP-43 pathology is unclassifiable")
level_vec[idx] <- "Unclassifiable"
levels(seurat_obj$`LATE-NC stage`) <- level_vec

head(seurat_obj@meta.data)

var_vec <- c("Age at death", "Cognitive status",
             "Braak stage", "Thal phase",
             "CERAD score", "APOE4 status",
             "LATE-NC stage")

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup2/Writeup2_microglia-pvm_umap-covariates.pdf",
    onefile = T, width = 5, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = "umap",
                           group.by = variable)
  print(plot1)
}

dev.off()

print("Done! :)")