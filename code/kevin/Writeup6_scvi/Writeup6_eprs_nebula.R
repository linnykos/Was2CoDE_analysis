rm(list=ls())

library(Seurat)

load("~/kzlinlab/data/jayadev_pu1-only_eprs/QCd_louvain_multi_clustered_ws_28pcs_ePRS_23s_uns_ADNChigh_20240503_anch.rdata")

# clear up workspace
ls_vec <- ls()
ls_vec <- setdiff(ls_vec, "ss_data_norm")
rm(list = ls_vec)
gc(TRUE)

tmp <- ss_data_norm$integrated_snn_res.0.1
tmp <- factor(paste0("c:", as.character(tmp)))
ss_data_norm$integrated_snn_res.0.1 <- tmp

# iterate over all cell clusters
cluster_names <- levels(ss_data_norm$integrated_snn_res.0.1)
nebula_res_list <- vector("list", length = length(cluster_names))
names(nebula_res_list) <- cluster_names

for(i in 1:length(cluster_names)){
  cluster_name <- cluster_names[i]
  print(Sys.time())
  print(paste("Working on cluster", cluster_name))
  
  ss_data_norm_subset <- subset(ss_data_norm, integrated_snn_res.0.1 == cluster_name)
  neb_data <- nebula::scToNeb(obj = ss_data_norm_subset,
                              assay = "RNA",
                              id = "Sample_ID",
                              pred = c("ePRS", "sex"),
                              offset = "nCount_RNA")
  df <- model.matrix(~ ePRS + sex,
                     data = neb_data$pred)
  if(is.unsorted(neb_data$id) == TRUE){
    neb_data <- nebula::group_cell(count = neb_data$count, 
                                   id = neb_data$id,
                                   pred = df,
                                   offset = neb_data$offset)
  } 
  
  nebula_res <- nebula::nebula(count = neb_data$count, 
                               id = neb_data$id,
                               pred = df,
                               offset = neb_data$offset,
                               model = "NBGMM",
                               verbose = TRUE)
  
  nebula_res_list[[i]] <- nebula_res
  
  save(nebula_res_list,
       file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_eprs_nebula_tmp.RData")
}

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Working from ~/kzlinlab/data/jayadev_pu1-only_eprs/QCd_louvain_multi_clustered_ws_28pcs_ePRS_23s_uns_ADNChigh_20240503_anch.rdata.",
              "Applying NEBULA.")

save(nebula_res, 
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_eprs_nebula.RData")

