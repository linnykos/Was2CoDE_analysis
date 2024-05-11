rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/jayadev_corbin_pbmc/QCd_louvain_multi_clustered_harmony_ws_24pcs_FAD_PBMC_20samples_20240321_harmony_by_sample.rdata")

# will only work if the assay is NOT in v5, see https://satijalab.org/seurat/articles/seurat5_essential_commands
class(ss_data_norm[["RNA"]])

# removing problmeatic things, to avoid the following error
# Error in slot(object = object, name = x) : 
# no slot of name "median_umi" for this object of class "SCTModel"
ss_data_norm[["SCT"]]@SCTModel.list <- list()

SeuratDisk::SaveH5Seurat(ss_data_norm,
                         filename = "~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_corbin_pbmc.h5Seurat",
                         overwrite = TRUE)
SeuratDisk::Convert("~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_corbin_pbmc.h5Seurat",
                    dest = "h5ad",
                    overwrite = TRUE)