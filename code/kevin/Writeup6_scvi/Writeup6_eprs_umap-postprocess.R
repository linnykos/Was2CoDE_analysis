rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)

load("~/kzlinlab/data/jayadev_pu1-only_eprs/QCd_norm_anchored_rpca__ws_anch_mg_subset.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

umap_csv <- read.csv("~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_eprs_mg_scvi_umap.csv")
rownames(umap_csv) <- umap_csv[,1]
umap_csv <- umap_csv[,-1]
colnames(umap_csv) <- paste0("scVI_UMAP_", 1:2)
umap_csv <- umap_csv[Seurat::Cells(ss_data_norm),]
umap_csv <- as.matrix(umap_csv)

ss_data_norm[["scvi_umap"]] <- Seurat::CreateDimReducObject(umap_csv,
                                                            assay = "integrated")
plot1 <- Seurat::DimPlot(ss_data_norm, 
                         reduction = "scvi_umap",
                         group.by = "Sample_ID",
                         raster = TRUE)
ggplot2::ggsave(plot1,
                filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_eprs_mg_scvi_Sample_ID_R.png",
                width = 8, height = 8)



