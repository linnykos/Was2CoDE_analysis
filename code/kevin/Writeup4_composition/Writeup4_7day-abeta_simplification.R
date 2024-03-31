rm(list=ls())
library(Seurat)
library(radEmu)

load("~/kzlinlab/data/ipsc_aqr7-2023/QCd_louvain_multi_clustered_harmony_ws_40pcs_20240205_SCTonly_7day_abeta_comp_run2.rdata")

keep_var <- "ss_data_norm"
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% keep_var]

rm(list = ls_vec)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("The data from ~/kzlinlab/data/ipsc_aqr7-2023/QCd_louvain_multi_clustered_harmony_ws_40pcs_20240205_SCTonly_7day_abeta_comp_run2.rdata,",
              "removing everything else, including ss_data_norm_unmerged (which is the",
              "data before merging)")

save(ss_data_norm, 
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_7day-abeta_simplified.RData")
