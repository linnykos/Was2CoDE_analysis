rm(list=ls())
library(Seurat)
library(fastTopics)
library(scCustomize)

load("~/kzlinlab/data/ipsc_aqr7-2023/QCd_louvain_multi_clustered_harmony_ws_40pcs_20240205_SCTonly_7day_abeta_comp_run2.rdata")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_7day-abeta_fasttopics.RData")
set.seed(10)

topic_mat <- matrix(NA, nrow = ncol(ss_data_norm), ncol = ncol(topic_res$L))
rownames(topic_mat) <- colnames(ss_data_norm)
topic_res$L <- topic_res$L[rownames(topic_res$L) %in% colnames(ss_data_norm),]

for(i in 1:nrow(topic_res$L)){
  if(i %% floor(nrow(topic_res$L)/10) == 0) cat('*')
  
  topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
}
colnames(topic_mat) <- paste0("fastTopic_", 1:ncol(topic_mat))
colnames(topic_res$F) <- paste0("fastTopic_", 1:ncol(topic_mat))

ss_data_norm[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                            loadings =  topic_res$F,
                                                            assay = "RNA",
                                                            key =  paste0("fastTopic_"))

######

num_per_plot <- 4

for(kk in 1:ceiling(ncol(topic_mat)/num_per_plot)){
  values <- ((kk-1)*num_per_plot+1):(min(kk*num_per_plot, ncol(topic_mat)))
  cat(c("Working on plots:", paste0(values, sep = ", "),"\n"))
  
  plot_list <- lapply(values, function(i){
    p1 <- scCustomize::FeaturePlot_scCustom(ss_data_norm, 
                                            colors_use = list("red", "lightgray", "blue"),
                                            na_color = "bisque",
                                            reduction = "umap", 
                                            features = paste0("fastTopic_", i))
    p1 <- p1 + ggplot2::ggtitle(paste0("UMAP, visualizing fastTopic ", i))
    
    p1
  })
  
  plot_all <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
  
  ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_7day-abeta_fasttopic",
                                    min(values), "-", max(values), ".png"),
                  plot_all, device = "png", 
                  width = 10, 
                  height = 5*round(length(plot_list)/2), 
                  units = "in")
}



