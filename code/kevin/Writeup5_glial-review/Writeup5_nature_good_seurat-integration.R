rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL

seurat_obj <- ss_data_norm
ls_vec <- ls(); ls_vec <- setdiff(ls_vec, seurat_obj)
rm(ls = ls_vec); gc(TRUE)

seurat_obj$Pt_ID <- paste0("D:", as.character(seurat_obj$Pt_ID))
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[which(seurat_obj$Pt_ID %in% c("D:1", "D:2", "D:3", "D:9", 
                                       "D:21", "D:22", "D:4", "D:5", 
                                       "D:6", "D:11"))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# # we will kick out some donors in this experiment
# donor_vec <- unique(metadf$Pt_ID)
# donor_df <- t(sapply(donor_vec, function(donor){
#   idx <- which(metadf$Pt_ID == donor)
#   as.character(metadf[idx[1], c("Pt_ID", "SeqBatch", "Study_Designation", "Sex")])
# }))
# donor_df <- as.data.frame(donor_df)
# colnames(donor_df) <- c("Pt_ID", "SeqBatch", "Study_Designation", "Sex")
# table(donor_df$SeqBatch, donor_df$Study_Designation)

# Commented out procedure (Since I've decided we should stick to as-close-as-possible to the original workflow)
# https://satijalab.org/seurat/articles/integration_introduction.html
# ss_data_norm[["RNA"]] <- split(ss_data_norm[["RNA"]], f = ss_data_norm$Pt_ID)
# set.seed(10)
# ss_data_norm <- Seurat::IntegrateLayers(object = ss_data_norm, 
#                                         method = Seurat::CCAIntegration, 
#                                         orig.reduction = "pca", 
#                                         new.reduction = "integrated.cca",
#                                         verbose = TRUE)
# 
# # re-join layers after integration
# ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

var_subset <- "SeqBatch" # originally, this was "Pt_ID"
print(paste("Subsetting dataset with", var_subset, sep = " "))
seurat_obj <- Seurat::SplitObject(seurat_obj, split.by = var_subset)
print(paste("Split Seurat object now has",
            length(unique(names(seurat_obj))), "groups.", sep = " "))

unwanted <- c("percent.mito")
print(paste("Running SCT on Seurat object and regressing out",
            unwanted, "variable.", sep = " "))
nfeatures <- 2000 # changing from the original 5000
normalize <- function (data) {
  data <- Seurat::SCTransform(data, 
                              variable.features.n = nfeatures,
                              do.scale = FALSE,
                              do.center = FALSE,
                              conserve.memory = TRUE,
                              return.only.var.genes = FALSE,
                              vars.to.regress = unwanted,
                              verbose = FALSE)
}
ss_data_norm <- lapply(1:length(seurat_obj), function(i){
  set.seed(10)
  print(paste0("Working on fold number ", i, " out of ", length(seurat_obj)))
  return(normalize(seurat_obj[[i]]))
})

rm(ls = "seurat_obj"); gc(TRUE)

save(ss_data_norm, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_seurat-integration_tmp.RData")

# Find the integration features similar across all samples for the top
# nFeatures of genes.
print("Finding Ifeatures")
Ifeatures <- Seurat::SelectIntegrationFeatures(ss_data_norm, 
                                               nfeatures = nfeatures)

# Ensure we're using the normalized data to integrate.
print("Prepping integration")
ss_data_norm <- Seurat::PrepSCTIntegration(ss_data_norm, 
                                           anchor.features = Ifeatures)

save(ss_data_norm, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_seurat-integration_tmp.RData")

print("Finding integration anchors using no reference samples.")
ss_data_norm <- Seurat::FindIntegrationAnchors(ss_data_norm,
                                               normalization.method = "SCT",
                                               reference = NULL,
                                               anchor.features = Ifeatures,
                                               reduction = "cca",
                                               dims = 1:30)

save(ss_data_norm, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_seurat-integration_tmp.RData")

print("Performing the integration")
ss_data_norm <- Seurat::IntegrateData(ss_data_norm,
                                      normalization.method = "SCT",
                                      dims = 1:30)

save(ss_data_norm, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_seurat-integration_tmp.RData")

print("Performing PCA and UMAP")
set.seed(10)
Seurat::DefaultAssay(ss_data_norm) <- "integrated"
# ss_data_norm <- Seurat::FindVariableFeatures(ss_data_norm, nfeatures = nfeatures)
ss_data_norm <- Seurat::RunPCA(ss_data_norm, 
                               features = Seurat::VariableFeatures(ss_data_norm), 
                               npcs = 50)
set.seed(10)
ss_data_norm <- Seurat::RunUMAP(object = ss_data_norm, dims = 1:30)

save(ss_data_norm, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_seurat-integration.RData")





