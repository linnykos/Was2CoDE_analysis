rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(arrow)

out_folder <- "~/kzlinlab/projects/subject-de/out/kevin/Writeup11/"
load(paste0(out_folder, "Writeup11_rosmap-initial.RData"))
seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)

# Reading the Feather file
df <- arrow::read_feather(paste0(out_folder, "Writeup11_rosmap_scvi.feather"))
rowname_vec <- as.matrix(df[,ncol(df)])[,1]
df <- df[,-ncol(df)]
df <- as.matrix(df)
rownames(df) <- rowname_vec

Seurat::DefaultAssay(seurat_obj) <- "RNA"

Seurat::VariableFeatures(seurat_obj) <- colnames(df)
seurat_obj <- subset(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
SeuratObject::LayerData(object = seurat_obj, 
                        assay = "RNA", 
                        layer = "data") <- t(df)

# compute PCA and UMAP
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             features = Seurat::VariableFeatures(seurat_obj),
                             verbose = FALSE)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:30)

##############################

# fixing all the defficiencies
# update the RNA assay
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

# put factors
factor_vars_list <- list(batch = TRUE,
                         celltype = TRUE,
                         Pt_ID = TRUE,
                         ADdiag3types = c("nonAD", "earlyAD", "lateAD"),
                         brainRegion = TRUE,
                         subject = TRUE,
                         orig.ident = TRUE,
                         seurat_clusters = TRUE,
                         Study = TRUE,
                         race = TRUE,
                         spanish = TRUE,
                         apoe_genotype = TRUE,
                         sex = TRUE,
                         APOEe4_status = TRUE,
                         braaksc = TRUE,
                         ceradsc = TRUE,
                         cogdx = TRUE,
                         dcfdx_lv = TRUE,
                         ADpath = c("no", "yes"))
for(kk in 1:length(factor_vars_list)){
  variable <- names(factor_vars_list)[kk]
  col_idx <- which(colnames(seurat_obj@meta.data) == variable)
  vec <- as.character(seurat_obj@meta.data[,col_idx])
  if(all(factor_vars_list[[kk]] == TRUE)){
    level_vec <- sort(unique(vec))
  } else {
    level_vec <- factor_vars_list[[kk]]
  }
  seurat_obj@meta.data[,col_idx] <- factor(vec, levels = level_vec)
}

#######################

var_vec <- c(subject = FALSE,
             batch = TRUE,
             ADpath = TRUE,
             ADdiag3types = TRUE,
             seurat_clusters = TRUE,
             celltype = TRUE,
             Study = TRUE,
             race = TRUE,
             spanish = TRUE,
             apoe_genotype = TRUE,
             sex = TRUE,
             APOEe4_status = TRUE,
             braaksc = TRUE,
             ceradsc = TRUE,
             cogdx = TRUE,
             dcfdx_lv = TRUE)

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup11/Writeup11_rosmap_scvi_umap-covariates.pdf",
    onefile = T, width = 5, height = 5)

for(kk in 1:length(var_vec)){
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = "umap",
                           group.by = names(var_vec)[kk],
                           raster = TRUE)
  if(!var_vec[kk]) plot1 <- plot1 + Seurat::NoLegend()
  print(plot1)
}

dev.off()

#############

SeuratObject::LayerData(object = seurat_obj, 
                        assay = "RNA", 
                        layer = "scale.data") <- NULL

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- "After running scVI."
save(seurat_obj,
     date_of_run, session_info, note,
     file = paste0(out_folder, "Writeup11_rosmap_scVI-postprocessed.RData"))