rm(list=ls())

library(Seurat)
library(ideas)
library(IdeasCustom)
library(caret)
library(ggplot2)
library(data.table)
library(doRNG)
print(sessionInfo())

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

# the loaded dataset is called "ss_data_norm"

set.seed(10)

# https://satijalab.org/seurat/articles/essential_commands.html
count_matrix <- SeuratObject::LayerData(ss_data_norm,
                                       features = Seurat::VariableFeatures(ss_data_norm),
                                       layer = "data",
                                       assay = "RNA")

meta_cell <- ss_data_norm@meta.data
meta_cell$individual <- meta_cell$Pt_ID
meta_cell$cell_id <- row.names(meta_cell) 
meta_ind     <- data.frame("individual" = ss_data_norm$Pt_ID,
                           "Study_Designation" = ss_data_norm$Study_Designation,
                           "Study_Designation" = ss_data_norm$Study_Designation,
                           "Sex" = ss_data_norm$Sex,
                           "APOEe4_status" = ss_data_norm$APOEe4_status,
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

dmy_gen <- dummyVars("~ APOEe4_status", data = meta_ind)
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
                                  "Study_Designation" = ss_data_norm$Study_Designation,
                                  "Sex" = dmy_Sex,
                                  "APOEe4_status" = dmy_gen,
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

zz <- as.character(meta_ind[,"coded_Age"])
meta_ind[,"coded_Age"] <- as.numeric(zz)

for(j in 1:ncol(meta_ind)){
  if(!is.numeric(meta_ind[,j])) next()
  meta_ind[which(is.na(meta_ind[,j])),j] <- stats::median(meta_ind[,j], na.rm = T)
}

var2test      = "Study_Designation"
var2adjust    =  setdiff(colnames(meta_ind), c("individual", "Study_Designation"))
var2test_type = "binary" 
var_per_cell  =  "nCount_RNA" 

###########

save(meta_cell, meta_ind,
     var_per_cell, var2test, var2test_type, var2adjust,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_microglia_ideascustom.RData")

dist_list = IdeasCustom::ideas_dist_custom(count_input = count_matrix, 
                                           meta_cell = meta_cell, 
                                           meta_ind = meta_ind, 
                                           var_per_cell = var_per_cell, 
                                           var2test = var2test, 
                                           var2test_type = var2test_type,
                                           verbose = 3)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Tati's Was2 IDEAS analysis of the microglia data.",
              "This was done on the data in ~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData.")

save(meta_cell,
     meta_ind,
     dist_list,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_microglia_ideascustom.RData")

print("Done! :)")
