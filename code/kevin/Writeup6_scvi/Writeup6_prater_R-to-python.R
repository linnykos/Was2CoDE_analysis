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

# adjust the APOE meta-data column to avoid issues
col_idx <- which(colnames(ss_data_norm@meta.data) == "APOE")
colnames(ss_data_norm@meta.data)[col_idx] <- "genotype_APOE"
colnames(ss_data_norm@meta.data)

# change the CognitiveStatus levels to avoid any issues
ss_data_norm$CognitiveStatus <- factor(ss_data_norm$CognitiveStatus)
level_vec <- levels(ss_data_norm$CognitiveStatus)
level_vec[which(level_vec == "No dementia")] <- "No_dementia"
levels(ss_data_norm$CognitiveStatus) <- level_vec
table(ss_data_norm$CognitiveStatus)

# convert some variables into factors
variable_names <- c("Pt_ID", "Study_Designation", "Sex", "genotype_APOE", "Race")
for(variable in variable_names){
  ss_data_norm@meta.data[,variable] <- droplevels(factor(ss_data_norm@meta.data[,variable]))
}
summary(ss_data_norm@meta.data)

# fix NAs in race
vec <- as.character(ss_data_norm$Race)
vec[is.na(vec)] <- "NA"
ss_data_norm$Race <- factor(vec)
table(ss_data_norm$Race)

# fix coded_Age
vec <- as.character(ss_data_norm$coded_Age)
vec[vec == "90+"] <- "90"
ss_data_norm$coded_Age <- as.numeric(vec)
quantile(ss_data_norm$coded_Age)

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

# convert everything back into characters in preparation for the conversion
variable_names <- c("Pt_ID", "Study_Designation", "Sex", "genotype_APOE", "Race", 
                    "CognitiveStatus", "integrated_snn_res.0.3", "seurat_clusters", 
                    "orig.ident")
for(variable in variable_names){
  ss_data_norm@meta.data[,variable] <- as.character(ss_data_norm@meta.data[,variable])
}
summary(ss_data_norm@meta.data)

SeuratDisk::SaveH5Seurat(ss_data_norm, 
                         filename = "~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.h5Seurat")
SeuratDisk::Convert("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.h5Seurat", 
                    dest = "h5ad", 
                    misc = FALSE)

print("Done! :)")