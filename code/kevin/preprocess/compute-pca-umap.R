rm(list=ls())
library(Seurat)
set.seed(10)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

# see https://satijalab.org/seurat/articles/seurat5_essential_commands
options(Seurat.object.assay.version = "v5")

Seurat::DefaultAssay(ss_data_norm) <- "integrated"
ss_data_norm[["integrated"]] <- as(object = ss_data_norm[["integrated"]], Class = "Assay5")
ss_data_norm[["RNA"]] <- as(object = ss_data_norm[["RNA"]], Class = "Assay5")

# adjust the APOE meta-data column to avoid issues
col_idx <- which(colnames(ss_data_norm@meta.data) == "APOE")
colnames(ss_data_norm@meta.data)[col_idx] <- "genotype_APOE"

# change the CognitiveStatus levels to avoid any issues
ss_data_norm$CognitiveStatus <- factor(ss_data_norm$CognitiveStatus)
level_vec <- levels(ss_data_norm$CognitiveStatus)
level_vec[which(level_vec == "No dementia")] <- "No_dementia"
levels(ss_data_norm$CognitiveStatus) <- level_vec

# do a basic visualization
ss_data_norm <- Seurat::RunPCA(ss_data_norm, 
                               features = Seurat::VariableFeatures(ss_data_norm), 
                               npcs = 50)
ss_data_norm <- Seurat::RunUMAP(object = ss_data_norm, 
                                dims = 1:30)

# Saving file
date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Original data from ~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

save(ss_data_norm, 
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData")

print("Done! :)")
