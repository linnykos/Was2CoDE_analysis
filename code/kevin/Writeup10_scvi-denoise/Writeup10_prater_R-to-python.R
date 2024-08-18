rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

###################### # Preprocessing

# remove unnecessary objects
Seurat::DefaultAssay(ss_data_norm) <- "RNA"
ss_data_norm[["integrated"]] <- NULL
SeuratObject::LayerData(object = ss_data_norm, 
                        assay = "RNA", 
                        layer = "data") <- NULL

# adjust the APOE meta-data column to avoid issues
col_idx <- which(colnames(ss_data_norm@meta.data) == "APOE")
colnames(ss_data_norm@meta.data)[col_idx] <- "genotype_APOE"

# go through all the column names and remove all punctuations (aside from "_" and ".")
colnames_vec <- colnames(ss_data_norm@meta.data)
colnames_vec <- gsub(pattern = "[^[:alnum:]._ ]", replacement = "", x = colnames_vec)
colnames_vec <- gsub(pattern = " ", replacement = "", x = colnames_vec)
colnames(ss_data_norm@meta.data) <- colnames_vec

# go through each column. If it's a numeric, ignore. If it's a factor or character, remove all punctuations (aside from "_" and ".")
for(j in 1:ncol(ss_data_norm@meta.data)){
  tmp <- ss_data_norm@meta.data[,j]
  if(is.numeric(tmp)) next()
  if(is.factor(tmp)) tmp <- as.character(tmp)
  tmp <- gsub(pattern = "[^[:alnum:]._ ]", replacement = "", x = tmp)
  tmp <- gsub(pattern = " ", replacement = "", x = tmp)
  ss_data_norm@meta.data[,j] <- tmp
}

# add letters to certain clusterings
letters_vec <- c(SeqBatch = "b", 
                 seurat_clusters = "c", 
                 integrated_snn_res.0.3 = "c", 
                 Pt_ID = "pt", 
                 orig.ident ="pt",
                 genotype_APOE = "apoe")
for(k in 1:length(letters_vec)){
  variable <- names(letters_vec)[k]
  letter <- letters_vec[k]
  ss_data_norm@meta.data[,variable] <- paste0(letter, ss_data_norm@meta.data[,variable])
}

# need to fill in NAs in PMI,Race
# from: https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-023-00424-y/MediaObjects/43587_2023_424_MOESM1_ESM.pdf
pmi_missing_values <- c(pt1 = 5.08,
                        pt10 = 3.33,
                        pt3 = 5.08,
                        pt15 = 6.97)
for(i in 1:length(pmi_missing_values)){
  subj_id <- names(pmi_missing_values)[i]
  idx <- which(ss_data_norm$Pt_ID == subj_id)
  ss_data_norm$PMI[idx] <- pmi_missing_values[i]
}

# based on empirical matching analysis
race_missing_values <- c(pt5 = "White")
for(i in 1:length(race_missing_values)){
  subj_id <- names(race_missing_values)[i]
  idx <- which(ss_data_norm$Pt_ID == subj_id)
  ss_data_norm$Race[idx] <- race_missing_values[i]
}

# make sure there is age is cleaned
tmp <- as.character(ss_data_norm$coded_Age)
tmp[which(tmp == "90+")] <- "90"
ss_data_norm$coded_Age <- as.numeric(tmp)

# make sure all the numerical variables are numerics
numerical_vars <- c("coded_Age", "PMI")
for(variable in numerical_vars){
  ss_data_norm@meta.data[,variable] <- as.numeric(as.character(ss_data_norm@meta.data[,variable]))
}

# add an APOEe4_status column
tmp <- ss_data_norm$genotype_APOE
vec <- rep(0, length(tmp))
vec[which(tmp %in% c("apoe34", "apoe44"))] <- 1
ss_data_norm$APOEe4_status <- vec

# checking
categorical_vars <- c("Sex", "SeqBatch", "Race", "Study_Designation", "APOEe4_status", "Pt_ID")
numerical_vars <- c("PMI","coded_Age")

zz <- ss_data_norm@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(sapply(zz[,numerical_vars], is.numeric))
stopifnot(sapply(zz[,categorical_vars], is.numeric))
stopifnot(!any(is.na(zz)))
stopifnot(!any(zz == "NA"))
stopifnot(all(!sapply(1:ncol(zz), is.factor)))
summary(zz)
head(zz)

###################### # Saving as R

note <- "Doing some simple cleaning of Katie's data for better usage. Changed all factors into characters."
save(ss_data_norm,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_cleaned.RData")

###################### # Saving as python

SeuratDisk::SaveH5Seurat(ss_data_norm, 
                         filename = "~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_cleaned.h5Seurat")
SeuratDisk::Convert("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_cleaned.h5Seurat", 
                    dest = "h5ad", 
                    misc = FALSE)