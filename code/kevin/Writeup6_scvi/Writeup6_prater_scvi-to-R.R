rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(arrow)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

# Reading the Feather file
df <- arrow::read_feather('/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-noCovariates-denoised.feather')
rowname_vec <- as.matrix(df[,ncol(df)])[,1]
df <- df[,-ncol(df)]
df <- as.matrix(df)
rownames(df) <- rowname_vec

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL

# remove "scale.data" slot
SeuratObject::LayerData(object = ss_data_norm, 
                        assay = "RNA", 
                        layer = "scale.data") <- NULL

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

####################
# basic preprocessing

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["RNA"]] <- as(object = ss_data_norm[["RNA"]], Class = "Assay5")

# adjust the APOE meta-data column to avoid issues
col_idx <- which(colnames(ss_data_norm@meta.data) == "APOE")
colnames(ss_data_norm@meta.data)[col_idx] <- "genotype_APOE"

# change the CognitiveStatus levels to avoid any issues
ss_data_norm$CognitiveStatus <- factor(ss_data_norm$CognitiveStatus)
level_vec <- levels(ss_data_norm$CognitiveStatus)
level_vec[which(level_vec == "No dementia")] <- "No_dementia"
levels(ss_data_norm$CognitiveStatus) <- level_vec

# add letters to certain clusterings
variable_names <- c("SeqBatch", 
                    "seurat_clusters", 
                    "integrated_snn_res.0.3", 
                    "Pt_ID", 
                    "orig.ident")
letters_vec <- c("b", "c", "c", "pt", "pt")
for(k in 1:length(variable_names)){
  variable <- variable_names[k]
  letter <- letters_vec[k]
  ss_data_norm@meta.data[,variable] <- factor(paste0(letter, ss_data_norm@meta.data[,variable]))
}

# make things into their appropriate type
categorical_vars <- c("Sex", 
                      "Race", 
                      "CognitiveStatus", 
                      "genotype_APOE", 
                      "Study_Designation")
for(variable in categorical_vars){
  ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
}
summary(ss_data_norm@meta.data)

#######################

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- "This is Katie's data but putting the scVI-normalized data (only 2000 genes) inside under 'data'."

save(session_info, date_of_run, note,
     ss_data_norm, 
     file = "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData")

#######################

var_vec <- c("seurat_clusters", "Pt_ID", "CognitiveStatus", 
             "Study_Designation", "Sex", "genotype_APOE",
             "PMI", "BrainPh", "Race", "FreshBrainWeight",
             "NIA_AA", "ThalPhase", 
             "BraakStage", "CERAD", "LATEScore", "coded_Age",
             "SeqBatch")

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_prater_scvi-seurat_umap-covariates.pdf",
    onefile = T, width = 5, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(ss_data_norm, 
                           reduction = "umap",
                           group.by = variable,
                           raster = TRUE)
  print(plot1)
}

dev.off()

