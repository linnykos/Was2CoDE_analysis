rm(list=ls())
library(Seurat)

file_folder <- "~/kzlinlab/data/jayadev_corbin_pbmc/"
file_vec <- c("QCd_louvain_multi_clustered_harmony_ws_24pcs_FAD_PBMC_20samples_20240321_harmony_by_sample.rdata",
              "QCd_louvain_multi_clustered_harmony_ws_33pcs_monos_clust39101112.rdata",
              "QCd_louvain_multi_clustered_harmony_ws_32pcs_nktcells_clust012467814.rdata")

file_name <- c("24pcs_FAD_PBMC_20samples_20240321",
               "33pcs_monos",
               "32pcs_nktcells")
seurat_clustering <- c("0.25", "0.3", "0.25")

out_csv_folder <- "~/kzlinlab/projects/subject-de/out/kevin/Writeup4/"
out_plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup4/"

for(kk in 2:length(file_vec)){
  file <- file_vec[kk]
  print(paste0("Working on file: ", file))
  
  load(paste0(file_folder, file_vec[kk]))
  print("Done loading")
  
  metadata <- ss_data_norm@meta.data
  ss_data_norm$Pt_ID_mutation <- sapply(1:nrow(metadata), function(i){
    paste0(metadata[i,c("mutation_status", "Pt_ID")], collapse = "_")
  })
  
  Y <- table(ss_data_norm$Pt_ID_mutation, ss_data_norm@meta.data[,paste0("SCT_snn_res.", seurat_clustering[kk])])
  colnames(Y) <- paste0("cluster_", 0:(ncol(Y)-1))
  rowname_vec <- rownames(Y)
  Y <- as.data.frame(apply(as.matrix.noquote(Y),2,as.numeric))
  Y <- cbind(rowname_vec, Y)
  colnames(Y)[1] <- "Pt_ID_mutation"
  
  write.csv(Y, 
            row.names = FALSE,
            file = paste0(out_csv_folder, file_name[kk], "_composition.csv"))
  
  # other covariates
  other_covariate_vec <- c("sex", "AGE", "apoe_4_dose")
  donor_vec <- Y[,1]
  donor_metadf <- t(sapply(1:length(donor_vec), function(i){
    idx <- which(ss_data_norm$Pt_ID_mutation == donor_vec[i])[1]
    vec <- sapply(other_covariate_vec, function(variable){
      metadata[idx,variable]
    })
    return(as.character(c(donor_vec[i], vec)))
  }))
  donor_metadf <- as.data.frame(donor_metadf)
  colnames(donor_metadf) <- c("Pt_ID_mutation", other_covariate_vec)
  rownames(donor_metadf) <- NULL
  donor_metadf$AGE <- as.numeric(donor_metadf$AGE)
  donor_metadf$apoe_4_dose <- as.numeric(donor_metadf$apoe_4_dose)
  
  write.csv(donor_metadf, 
            row.names = FALSE,
            file = paste0(out_csv_folder, file_name[kk], "_donor-metadata.csv"))
  
  plot1 <- Seurat::DimPlot(ss_data_norm,
                           group.by = paste0("SCT_snn_res.", seurat_clustering[kk]))
  ggplot2::ggsave(filename = paste0(out_plot_folder, file_name[kk], "_umap.png"),
                  plot1, device = "png", width = 8, height = 5, units = "in",
                  dpi = 300)
  
  ##########
  
  plot1 <- Seurat::DimPlot(ss_data_norm,
                           group.by = "mutation_status")
  plot2 <- Seurat::DimPlot(ss_data_norm,
                           group.by = "sex")
  plot3 <- Seurat::DimPlot(ss_data_norm,
                           group.by = "apoe_4_dose")
  
  plot_all <- cowplot::plot_grid(plot1, plot2, plot3, ncol = 3)
  ggplot2::ggsave(filename = paste0(out_plot_folder, file_name[kk], "_umap-covariates.png"),
                  plot_all, device = "png", width = 15, height = 5, units = "in",
                  dpi = 300)
}

print("Done! :)")