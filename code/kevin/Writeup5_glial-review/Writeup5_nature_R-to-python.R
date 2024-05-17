rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL
ss_data_norm$Pt_ID <- paste0("D:", as.character(ss_data_norm$Pt_ID))

table(ss_data_norm$Pt_ID)

SeuratObject::LayerData(object = ss_data_norm, 
                        assay = "RNA", 
                        layer = "data") <- NULL

tmp <- SeuratObject::LayerData(object = ss_data_norm, 
                               assay = "RNA", 
                               layer = "counts") 
head(tmp@x)

tmp <- SeuratObject::LayerData(object = ss_data_norm, 
                               assay = "RNA", 
                               layer = "data") 
head(tmp@x)

SeuratDisk::SaveH5Seurat(ss_data_norm, 
                         filename = "~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.h5Seurat")
SeuratDisk::Convert("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.h5Seurat", 
                    dest = "h5ad", 
                    misc = FALSE)