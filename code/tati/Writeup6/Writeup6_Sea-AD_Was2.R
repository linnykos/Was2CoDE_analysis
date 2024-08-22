rm(list=ls())

library(Seurat)
library(ideas)
library(IdeasCustom)
library(caret)
library(ggplot2)
library(data.table)
library(doRNG)
print(sessionInfo())

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData") 
# the loaded dataset is called "seurat_obj"

set.seed(10)

# https://satijalab.org/seurat/articles/essential_commands.html
count_matrix <- SeuratObject::LayerData(seurat_obj,
                                        features = Seurat::VariableFeatures(seurat_obj),
                                        layer = "data",
                                        assay = "RNA")

meta_cell <- seurat_obj@meta.data
meta_cell$individual <- meta_cell$donor_id
meta_cell$cell_id <- row.names(meta_cell) 
meta_ind     <- data.frame("individual" = seurat_obj$donor_id,
                           "ADNC" = seurat_obj$ADNC,
                           "sex" = seurat_obj$sex,
                           "APOE4status" = seurat_obj$APOE4status,
                           "PMI" = seurat_obj$PMI,
                           "BrainPh" = seurat_obj$BrainPh,
                           "Race" = seurat_obj$Race ,
                           "FreshBrainWeight" = seurat_obj$FreshBrainWeight,
                           "NIA_AA" = seurat_obj$NIA_AA,
                           "ThalPhase" = seurat_obj$ThalPhase,
                           "BraakStage" = seurat_obj$BraakStage,
                           "CERAD" = seurat_obj$CERAD,
                           "LATEScore" = seurat_obj$LATEScore,
                           "SeqBatch" = seurat_obj$SeqBatch,
                           "coded_Age" = seurat_obj$coded_Age,
                           row.names=NULL
)
dmy_SD<- dummyVars("~ADNC", data = meta_ind)
dmy_SD <- data.frame(predict(dmy_SD, newdata = meta_ind))
dmy_SD<- dmy_SD[,-which.max(colSums(dmy_SD)),drop=FALSE]
dmy_sex <- dummyVars("~sex", data = meta_ind)
dmy_sex <- data.frame(predict(dmy_sex, newdata = meta_ind))
dmy_sex<- dmy_sex[,-which.max(colSums(dmy_sex)),drop=FALSE]
dmy_Race <- dummyVars("~Race", data = meta_ind)
dmy_Race <- data.frame(predict(dmy_Race, newdata = meta_ind))

for(j in 1:ncol(dmy_Race)){
  dmy_Race[which(is.na(dmy_Race[,j])),j] <- 0
}
dmy_Race<- dmy_Race[,-which.max(colSums(dmy_Race)),drop=FALSE]

dmy_gen <- dummyVars("~ APOE4status", data = meta_ind)
dmy_gen <- data.frame(predict(dmy_gen, newdata = meta_ind))
dmy_gen<- dmy_gen[,-which.max(colSums(dmy_gen)),drop=FALSE]

meta_ind$SeqBatch <- factor(meta_ind$SeqBatch)
dmy_SB<- dummyVars("~SeqBatch", data = meta_ind)
dmy_SB <- data.frame(predict(dmy_SB, newdata = meta_ind))

for(j in 1:ncol(dmy_SB)){
  dmy_SB[which(is.na(dmy_SB[,j])),j] <- 0
}
dmy_SB<- dmy_SB[,-which.max(colSums(dmy_SB)),drop=FALSE]

meta_ind     <- unique(data.frame("individual" = seurat_obj$donor_id,
                                  "ADNC" = dmy_SD,
                                  "CognitiveStatus" = seurat_obj$CognitiveStatus,
                                  "sex" = dmy_sex,
                                  "APOE4status" = dmy_gen,
                                  "PMI" = seurat_obj$PMI,
                                  "BrainPh" = seurat_obj$BrainPh,
                                  "Race" = dmy_Race,
                                  "FreshBrainWeight" = seurat_obj$FreshBrainWeight,
                                  "NIA_AA" = seurat_obj$NIA_AA,
                                  "ThalPhase" = seurat_obj$ThalPhase,
                                  "BraakStage" = seurat_obj$BraakStage,
                                  "CERAD" = seurat_obj$CERAD,
                                  "LATEScore" = seurat_obj$LATEScore,
                                  "SeqBatch" = seurat_obj$SeqBatch,
                                  "coded_Age" = seurat_obj$coded_Age,
                                  row.names=NULL
))

zz <- as.character(meta_ind[,"coded_Age"])
zz[zz == "90+"] <- "90"
meta_ind[,"coded_Age"] <- as.numeric(zz)

for(j in 1:ncol(meta_ind)){
  if(!is.numeric(meta_ind[,j])) next()
  meta_ind[which(is.na(meta_ind[,j])),j] <- stats::median(meta_ind[,j], na.rm = T)
}

var2test      = "CognitiveStatus"
var2adjust    =  setdiff(colnames(meta_ind), c("individual", "CognitiveStatus"))
var2test_type = "binary" 
var_per_cell  =  "nCount_RNA" 

###########

save(meta_cell, meta_ind,
     var_per_cell, var2test, var2test_type, var2adjust,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_ideascustom.RData")

# load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup_microglia_ideascustom.RData")
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
              "This was done on the data in ~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData.")

save(meta_cell,
     meta_ind,
     dist_list,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup6_SEA-ADmicroglia_ideascustom.RData")

print("Done! :)")