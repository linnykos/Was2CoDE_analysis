rm(list=ls())
library(Seurat)
library(fastTopics)

load("~/kzlinlab/data/ipsc_aqr7-2023/QCd_louvain_multi_clustered_harmony_ws_40pcs_20240205_SCTonly_7day_abeta_comp_run2.rdata")
set.seed(10)

mat <- Seurat::GetAssayData(object = ss_data_norm, 
                            assay = "RNA", 
                            layer = "counts")
mat <- mat[Seurat::VariableFeatures(ss_data_norm),]
mat <- Matrix::t(mat)

K <- 10
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(topic_res, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_7day-abeta_fasttopics.RData")

print("Done! :)")

