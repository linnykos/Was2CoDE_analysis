rm(list=ls())
library(Seurat)
set.seed(10)

load("~/lab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

# see https://satijalab.org/seurat/articles/seurat5_essential_commands
options(Seurat.object.assay.version = "v5")

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL
ss_data_norm[["RNA"]] <- as(object = ss_data_norm[["RNA"]], Class = "Assay5")

# adjust the APOE meta-data column to avoid issues
col_idx <- which(colnames(ss_data_norm@meta.data) == "APOE")
colnames(ss_data_norm@meta.data)[col_idx] <- "genotype_APOE"

# zz <- ss_data_norm[["RNA"]]
# class(zz)
# names(zz@layers)
# zz@layers$counts[1:5,1:5]
# 
# # see https://rdrr.io/github/mojaveazure/seurat-object/src/R/logmap.R
# class(zz@cells)
# head(as.matrix(zz@cells))
# 
# class(zz@features)
# head(as.matrix(zz@features))
# 
# head(Seurat::Cells(ss_data_norm))
# head(colnames(ss_data_norm))

ss_data_norm <- Seurat::NormalizeData(object = ss_data_norm)
ss_data_norm <- Seurat::FindVariableFeatures(object = ss_data_norm)
ss_data_norm <- Seurat::ScaleData(object = ss_data_norm)
ss_data_norm <- Seurat::RunPCA(object = ss_data_norm)
ss_data_norm <- Seurat::RunUMAP(object = ss_data_norm, dims = 1:30)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Original data from ~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata.",
              "Trying with no integration at all.")

save(ss_data_norm, 
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/kevin/preprocess/naive-preprocess.RData")

print("Done! :)")

