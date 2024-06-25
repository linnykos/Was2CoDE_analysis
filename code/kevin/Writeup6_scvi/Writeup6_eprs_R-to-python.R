rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)

load("~/kzlinlab/data/jayadev_pu1-only_eprs/QCd_norm_anchored_rpca__ws_anch_mg_subset.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL
table(ss_data_norm$Sample_ID)

# append letters to all the clusterings
variable_names <- colnames(ss_data_norm@meta.data)[grep("integrated_snn", colnames(ss_data_norm@meta.data))]
variable_names <- c(variable_names, "seurat_clusters")
for(variable in variable_names){
  ss_data_norm@meta.data[,variable] <- paste0("c", as.character(ss_data_norm@meta.data[,variable]))
}

# append leters to all batches
variable_names <- colnames(ss_data_norm@meta.data)[grep("batch", colnames(ss_data_norm@meta.data))]
for(variable in variable_names){
  ss_data_norm@meta.data[,variable] <- paste0("b", as.character(ss_data_norm@meta.data[,variable]))
}
table(ss_data_norm$seq_batch, ss_data_norm$GEM_batch)
table(ss_data_norm$seq_batch, ss_data_norm$Library_batch)
summary(ss_data_norm@meta.data)
head(ss_data_norm@meta.data)

# remove "data" slot
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
                         filename = "~/kzlinlab/data/jayadev_pu1-only_eprs/QCd_norm_anchored_rpca__ws_anch_mg_subset.h5Seurat")
SeuratDisk::Convert("~/kzlinlab/data/jayadev_pu1-only_eprs/QCd_norm_anchored_rpca__ws_anch_mg_subset.h5Seurat", 
                    dest = "h5ad", 
                    misc = FALSE)

print("Done! :)")