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

set.seed(10)

count_matrix <- SeuratObject::LayerData(seurat_obj,
                                        features = Seurat::VariableFeatures(seurat_obj),
                                        layer = "data",
                                        assay = "RNA")

meta_cell <- seurat_obj@meta.data
meta_cell$individual <- meta_cell$donor_id
meta_cell$cell_id <- row.names(meta_cell) 

meta_ind     <- unique(data.frame("individual" = seurat_obj$donor_id,
                                  row.names=NULL
))


for(j in 1:ncol(meta_ind)){
  if(!is.numeric(meta_ind[,j])) next()
  meta_ind[which(is.na(meta_ind[,j])),j] <- stats::median(meta_ind[,j], na.rm = T)
}

var2test      = "ADNC"
var2adjust    =  setdiff(colnames(meta_ind), c("individual", "ADNC"))
var2test_type = "binary" 
var_per_cell  =  "nCount_RNA" 

###########

save(meta_cell, meta_ind,
     var_per_cell, var2test, var2test_type, var2adjust,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_ideascustom.RData")

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
              "This was done on the data in ~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData.")

save(meta_cell,
     meta_ind,
     dist_list,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_microglia_ideascustom.RData")

print("Done! :)")