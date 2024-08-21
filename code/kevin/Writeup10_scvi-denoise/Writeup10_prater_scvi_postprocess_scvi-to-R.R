rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(arrow)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_cleaned.RData")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

# Reading the Feather file
df <- arrow::read_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scvi.feather")
rowname_vec <- as.matrix(df[,ncol(df)])[,1]
df <- df[,-ncol(df)]
df <- as.matrix(df)
rownames(df) <- rowname_vec

# df2 <- scale(df); cor_mat <- crossprod(df2)/nrow(df2); diag(cor_mat) <- NA

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL

Seurat::VariableFeatures(ss_data_norm) <- colnames(df)
ss_data_norm <- subset(ss_data_norm, features = Seurat::VariableFeatures(ss_data_norm))
SeuratObject::LayerData(object = ss_data_norm, 
                        assay = "RNA", 
                        layer = "data") <- t(df)

# compute PCA and UMAP
ss_data_norm <- Seurat::ScaleData(ss_data_norm)
ss_data_norm <- Seurat::RunPCA(ss_data_norm, 
                               features = Seurat::VariableFeatures(ss_data_norm),
                               verbose = FALSE)
set.seed(10)
ss_data_norm <- Seurat::RunUMAP(ss_data_norm, dims = 1:30)

##############################

# fixing all the defficiencies
# update the RNA assay
ss_data_norm[["RNA"]] <- as(object = ss_data_norm[["RNA"]], Class = "Assay5")

# change 1 in "Braak" to "BraakI", 1 in "ThalPhase" is "Thal1"
# LATEScore is "LateStage1", etc..
# CERAD is "Absent", "Sparse", "Frequent", "Moderate"
cov_list <- list(
  BraakStage = c(Braak0 = 0, BraakI = 1, BraakII = 2, BraakIII = 3, BraakIV = 4, BraakV = 5, BraakVI = 6),
  CERAD = c(Absent = 0, Sparse = 1, Frequent = 2, Moderate = 3),
  LATEScore = c(LateStage0 = 0, LateStage1 = 1, LateStage2 = 2),
  NIA_AA = c(NIA0 = 0, NIA1 = 1, NIA2 = 2, NIA3 = 3), 
  ThalPhase = c(Thal0 = 0, Thal1 = 1, Thal2 = 2, Thal3 = 3, Thal4 = 4, Thal5 = 5)
)
for(kk in 1:length(cov_list)){
  col_idx <- which(colnames(ss_data_norm@meta.data) == names(cov_list)[kk])
  template <- cov_list[[kk]]
  vec <- as.numeric(ss_data_norm@meta.data[,col_idx])
  vec <- plyr::mapvalues(vec, from = template, to = names(template))
  ss_data_norm@meta.data[,col_idx] <- vec
}

# modification of APOEe4_status
vec <- ss_data_norm@meta.data[,"APOEe4_status"]
ss_data_norm@meta.data[,"APOEe4_status"] <- plyr::mapvalues(vec, from = c(0,1), to = c("N", "Y"))

# put factors
factor_vars_list <- list(BraakStage = TRUE,
                         CERAD = c("Absent", "Sparse", "Frequent", "Moderate"),
                         LATEScore = TRUE,
                         NIA_AA = TRUE,
                         ThalPhase = TRUE,
                         Sex = TRUE,
                         SeqBatch = TRUE,
                         Race = TRUE,
                         CognitiveStatus = c("Nodementia", "Dementia"),
                         Study_Designation = c("Ctrl", "AD"),
                         APOEe4_status = c("N", "Y"),
                         Pt_ID = names(sort(table(ss_data_norm$Pt_ID), decreasing = TRUE)),
                         seurat_clusters = TRUE,
                         integrated_snn_res.0.3 = TRUE,
                         genotype_APOE = TRUE,
                         orig.ident = TRUE)
for(kk in 1:length(factor_vars_list)){
  variable <- names(factor_vars_list)[kk]
  col_idx <- which(colnames(ss_data_norm@meta.data) == variable)
  vec <- as.character(ss_data_norm@meta.data[,col_idx])
  if(all(factor_vars_list[[kk]] == TRUE)){
    level_vec <- sort(unique(vec))
  } else {
    level_vec <- factor_vars_list[[kk]]
  }
  ss_data_norm@meta.data[,col_idx] <- factor(vec, levels = level_vec)
}

#######################

var_vec <- c("seurat_clusters", "Pt_ID", 
             "CognitiveStatus", "Study_Designation", 
             "Sex", "Race", "SeqBatch", "APOEe4_status",
             "PMI", "coded_Age",
             "NIA_AA", "ThalPhase", "CERAD", "LATEScore", "BraakStage")

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup10/Writeup10_prater_scvi_umap-covariates.pdf",
    onefile = T, width = 5, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(ss_data_norm, 
                           reduction = "umap",
                           group.by = variable,
                           raster = TRUE)
  print(plot1)
}

dev.off()


