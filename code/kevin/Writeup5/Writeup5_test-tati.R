# https://stackoverflow.com/questions/21171142/how-to-install-r-package-from-private-repo-using-devtools-install-github
# https://stackoverflow.com/questions/66065099/how-to-update-github-authentification-token-on-rstudio-to-match-the-new-policy
# usethis::use_git_config(user.name = "linnykos", 
#                         user.email = "kzlin@uw.edu")
# 
# credentials::set_github_pat()
# 
# devtools::install_github("TatiZhang/IdeasCustom",
#                          ref = "main")

#########

m(list=ls())

library(Seurat)
library(ideas)
library(IdeasCustom)
library(caret)
library(ggplot2)
library(data.table)
library(doRNG)

load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData") 
set.seed(10)

# https://satijalab.org/seurat/articles/essential_commands.html
count_matrix = SeuratObject::LayerData(ss_data_norm,
                                       features = Seurat::VariableFeatures(ss_data_norm),
                                       layer = "scale.data",
                                       assay = "integrated")
meta_cell    = ss_data_norm@meta.data
meta_cell$individual <- meta_cell$Pt_ID
meta_cell$cell_id <- row.names(meta_cell) 
meta_ind     <- data.frame("individual" = ss_data_norm$Pt_ID,
                           "Study_Designation" = ss_data_norm$Study_Designation,
                           "CognitiveStatus" = ss_data_norm$CognitiveStatus,
                           "Sex" = ss_data_norm$Sex,
                           "genotype_APOE" = ss_data_norm$genotype_APOE,
                           "PMI" = ss_data_norm$PMI,
                           "BrainPh" = ss_data_norm$BrainPh,
                           "Race" = ss_data_norm$Race ,
                           "FreshBrainWeight" = ss_data_norm$FreshBrainWeight,
                           "NIA_AA" = ss_data_norm$NIA_AA,
                           "ThalPhase" = ss_data_norm$ThalPhase,
                           "BraakStage" = ss_data_norm$BraakStage,
                           "CERAD" = ss_data_norm$CERAD,
                           "LATEScore" = ss_data_norm$LATEScore,
                           "SeqBatch" = ss_data_norm$SeqBatch,
                           "coded_Age" = ss_data_norm$coded_Age,
                           row.names=NULL
)
dmy_SD<- dummyVars("~Study_Designation", data = meta_ind)
dmy_SD <- data.frame(predict(dmy_SD, newdata = meta_ind))
dmy_SD<- dmy_SD[,-which.max(colSums(dmy_SD)),drop=FALSE]
dmy_Sex <- dummyVars("~Sex", data = meta_ind)
dmy_Sex <- data.frame(predict(dmy_Sex, newdata = meta_ind))
dmy_Sex<- dmy_Sex[,-which.max(colSums(dmy_Sex)),drop=FALSE]
dmy_Race <- dummyVars("~Race", data = meta_ind)
dmy_Race <- data.frame(predict(dmy_Race, newdata = meta_ind))

for(j in 1:ncol(dmy_Race)){
  dmy_Race[which(is.na(dmy_Race[,j])),j] <- 0
}
dmy_Race<- dmy_Race[,-which.max(colSums(dmy_Race)),drop=FALSE]

dmy_gen <- dummyVars("~ genotype_APOE", data = meta_ind)
dmy_gen <- data.frame(predict(dmy_gen, newdata = meta_ind))
dmy_gen<- dmy_gen[,-which.max(colSums(dmy_gen)),drop=FALSE]

meta_ind$SeqBatch <- factor(meta_ind$SeqBatch)
dmy_SB<- dummyVars("~SeqBatch", data = meta_ind)
dmy_SB <- data.frame(predict(dmy_SB, newdata = meta_ind))

for(j in 1:ncol(dmy_SB)){
  dmy_SB[which(is.na(dmy_SB[,j])),j] <- 0
}
dmy_SB<- dmy_SB[,-which.max(colSums(dmy_SB)),drop=FALSE]

meta_ind     <- unique(data.frame("individual" = ss_data_norm$Pt_ID,
                                  "Study_Designation" = dmy_SD,
                                  "CognitiveStatus" = ss_data_norm$CognitiveStatus,
                                  "Sex" = dmy_Sex,
                                  "genotype_APOE" = dmy_gen,
                                  "PMI" = ss_data_norm$PMI,
                                  "BrainPh" = ss_data_norm$BrainPh,
                                  "Race" = dmy_Race,
                                  "FreshBrainWeight" = ss_data_norm$FreshBrainWeight,
                                  "NIA_AA" = ss_data_norm$NIA_AA,
                                  "ThalPhase" = ss_data_norm$ThalPhase,
                                  "BraakStage" = ss_data_norm$BraakStage,
                                  "CERAD" = ss_data_norm$CERAD,
                                  "LATEScore" = ss_data_norm$LATEScore,
                                  "SeqBatch" = ss_data_norm$SeqBatch,
                                  "coded_Age" = ss_data_norm$coded_Age,
                                  row.names=NULL
))

class(meta_ind[,"coded_Age"])
levels(meta_ind[,"coded_Age"])
zz <- as.character(meta_ind[,"coded_Age"])
zz[zz == "90+"] <- "90"
meta_ind[,"coded_Age"] <- as.numeric(zz)

table(is.na(meta_ind))
for(j in 1:ncol(meta_ind)){
  if(!is.numeric(meta_ind[,j])) next()
  meta_ind[which(is.na(meta_ind[,j])),j] <- stats::median(meta_ind[,j], na.rm = T)
}

var2test      = "CognitiveStatus"
var2adjust    =  setdiff(colnames(meta_ind), c("individual", "CognitiveStatus"))
var2test_type = "binary" 
var_per_cell  =  "nCount_SCT" 

count_matrix_subset <- count_matrix[1:6,]
count_matrix_subset = as.matrix(count_matrix_subset)

dist_list = IdeasCustom::ideas_dist_custom(count_input = count_matrix, 
                                           meta_cell = meta_cell, 
                                           meta_ind = meta_ind, 
                                           var_per_cell = var_per_cell, 
                                           var2test = var2test, 
                                           var2test_type = var2test_type,
                                           verbose = 3)
