rm(list=ls())

library(Seurat)
library(ideas)
load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData") 
# the loaded dataset is called "ss_data_norm"
set.seed(10)

# https://satijalab.org/seurat/articles/essential_commands.html
count_matrix = SeuratObject::LayerData(ss_data_norm,
                                       features = Seurat::VariableFeatures(ss_data_norm),
                                       layer = "counts",
                                       assay = "RNA")
meta_cell    = ss_data_norm@meta.data
meta_ind     <- unique(data.frame("Individual" = ss_data_norm$Pt_ID, 
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
                 ))

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
var2adjust    =  c("Sex", "Race", "SeqBatch", "coded_Age") # [[KL: Thes are covariates you want to adjust for. It's not always obvious what to include here]]
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

count_matrix_subset <- count_matrix[1:10,]

dist1 = ideas_dist(count_matrix_subset, meta_cell, meta_ind, 
                   var_per_cell, var2test, var2test_type, 
                   d_metric = "Was", fit_method = "nb")

pval_ideas = permanova(dist1, meta_ind, var2test, var2adjust, 
                       var2test_type, n_perm=999, r.seed=903)


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
