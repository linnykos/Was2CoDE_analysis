rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/ipsc_aqr7-2023/QCd_louvain_multi_clustered_harmony_ws_40pcs_20240205_SCTonly_7day_abeta_comp_run2.rdata")

# will only work if the assay is NOT in v5, see https://satijalab.org/seurat/articles/seurat5_essential_commands
class(ss_data_norm[["RNA"]])

# removing problmeatic things, to avoid the following error
# Error in slot(object = object, name = x) : 
# no slot of name "median_umi" for this object of class "SCTModel"
ss_data_norm[["SCT"]]@SCTModel.list <- list()

SeuratDisk::SaveH5Seurat(ss_data_norm,
                         filename = "~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_7day-abeta.h5Seurat",
                         overwrite = TRUE)
SeuratDisk::Convert("~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_7day-abeta.h5Seurat",
                    dest = "h5ad",
                    overwrite = TRUE)