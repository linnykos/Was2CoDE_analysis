rm(list=ls())

library(Seurat)
library(ideas)
library(caret)
library(ggplot2)
load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData") 
# the loaded dataset is called "ss_data_norm"
set.seed(10)

# https://satijalab.org/seurat/articles/essential_commands.html
count_matrix = SeuratObject::LayerData(ss_data_norm,
                                       features = Seurat::VariableFeatures(ss_data_norm),
                                       layer = "counts",
                                       assay = "RNA")
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
str(meta_ind)
dmy_SD<- dummyVars("~Study_Designation", data = meta_ind)
dmy_SD <- data.frame(predict(dmy_SD, newdata = meta_ind))
head(dmy_SD)
dmy_SD<- dmy_SD[,-which.max(colSums(dmy_SD)),drop=FALSE]
str(meta_ind)
dmy_Sex <- dummyVars("~Sex", data = meta_ind)
dmy_Sex <- data.frame(predict(dmy_Sex, newdata = meta_ind))
dmy_Sex<- dmy_Sex[,-which.max(colSums(dmy_Sex)),drop=FALSE]
head(dmy_Sex)
dmy_Race <- dummyVars("~Race", data = meta_ind)
dmy_Race <- data.frame(predict(dmy_Race, newdata = meta_ind))
for(j in 1:ncol(dmy_Race)){
  dmy_Race[which(is.na(dmy_Race[,j])),j] <- 0
}
dmy_Race<- dmy_Race[,-which.max(colSums(dmy_Race)),drop=FALSE]
head(dmy_Race)
dmy_gen <- dummyVars("~ genotype_APOE", data = meta_ind)
dmy_gen <- data.frame(predict(dmy_gen, newdata = meta_ind))
dmy_gen<- dmy_gen[,-which.max(colSums(dmy_gen)),drop=FALSE]
head(dmy_gen)
meta_ind$SeqBatch <- factor(meta_ind$SeqBatch)
dmy_SB<- dummyVars("~SeqBatch", data = meta_ind)
dmy_SB <- data.frame(predict(dmy_SB, newdata = meta_ind))
for(j in 1:ncol(dmy_SB)){
  dmy_SB[which(is.na(dmy_SB[,j])),j] <- 0
}
dmy_SB<- dmy_SB[,-which.max(colSums(dmy_SB)),drop=FALSE]
head(dmy_SB)

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
dim(meta_ind)

# meta_ind     <- unique(data.frame("individual" = ss_data_norm$Pt_ID, 
#                                   "CognitiveStatus" = ss_data_norm$CognitiveStatus,
#                                   "coded_Age" = ss_data_norm$coded_Age,
#                                   row.names=NULL
# ))

class(meta_ind[,"coded_Age"])
levels(meta_ind[,"coded_Age"])
zz <- as.character(meta_ind[,"coded_Age"])
zz[zz == "90+"] <- "90"
meta_ind[,"coded_Age"] <- as.numeric(zz)
summary(meta_ind)

table(is.na(meta_ind))
for(j in 1:ncol(meta_ind)){
  if(!is.numeric(meta_ind[,j])) next()
  meta_ind[which(is.na(meta_ind[,j])),j] <- stats::median(meta_ind[,j], na.rm = T)
}

###########
# a quick side-demo on how factors-to-numerics can cause a lot of bugs

tmp <- factor(c("5","5","5","2","10","5","10"))
as.numeric(tmp) # gives you numbers according to the levels
as.numeric(as.character(tmp))

tmp <- factor(c("sadf", "egg", "egg", "sadf", "egg", "potato", "statistics", "sadf"))
as.numeric(tmp) # gives you numbers according to the levels
as.numeric(as.character(tmp))

  ###########

  # the main one to fill in. It should include the following: (one row per Pt_ID)
#   Study_Designation
# CognitiveStatus
# Sex
# genotype_APOE
# PMI
# BrainPh
# Race
# FreshBrainWeight
# NIA_AA
# ThalPhase
# BraakStage
# CERAD
# LATEScore
# SeqBatch
# coded_Age



var2test      = "CognitiveStatus"
var2adjust    =  setdiff(colnames(meta_ind), c("individual", "CognitiveStatus"))
# c("Sex", "Race", "SeqBatch", "coded_Age") # [[KL: Thes are covariates you want to adjust for. It's not always obvious what to include here]]
# var2adjust = "coded_Age"
var2test_type = "binary" # [[KL: In general, keep this as "binary"]]
var_per_cell  =  "nCount_SCT" # [[KL: This is the read depth, don't worry about this. I will tell you what to put here]]

###########
# 
# unique(ss_data_norm$Pt_ID)
# length(unique(ss_data_norm$Pt_ID))
# 
# dim(count_matrix)
# count_matrix[1:5,1:5]
# 
# str_to_factor_vec <- c("Study_Designation", "Sex", "genotype_APOE", "Race", "SeqBatch", "NIA_AA")
# for(variable in str_to_factor_vec){
#   ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
# }
# meta_cell    = ss_data_norm@meta.data
# head(meta_cell)
# summary(meta_cell)
# 
# table(ss_data_norm$Pt_ID,ss_data_norm$CognitiveStatus)
# table(ss_data_norm$Pt_ID,ss_data_norm$Sex)
# table(ss_data_norm$CognitiveStatus,ss_data_norm$Sex)

############

# count_matrix_subset <- count_matrix[1:10,]
# count_matrix_subset = as.matrix(count_matrix_subset)
# dist1 = ideas_dist(count_matrix_subset, meta_cell, meta_ind, 
#                    var_per_cell, var2test, var2test_type, 
#                    d_metric = "Was", fit_method = "kde")

count_matrix <- as.matrix(count_matrix)

dist1 = ideas_dist(count_matrix, meta_cell, meta_ind, 
                   var_per_cell, var2test, var2test_type, 
                   d_metric = "Was", fit_method = "kde")

save(dist1,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideas_tmp.RData")

pval_ideas = permanova(dist1, meta_ind, var2test, var2adjust, 
                       var2test_type, n_perm=999, r.seed=903)
head(pval_ideas)

##################

zz = meta_ind[,var2adjust]

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Basic IDEAS analysis of the microglia data.")

save(meta_cell,
     meta_ind,
     dist1,
     pval_ideas,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideas.RData")

print("Done! :)")
