rm(list=ls())
library(Seurat)
set.seed(10)

file_folder <- "~/kzlinlab/data/ipsc_aqr7-2023/"
file_vec <- c("QCd_louvain_multi_clustered_harmony_ws_40pcs_20240112_SCTonly_mgsubset_24hr_v_7d_abeta_comp.rdata",
              "QCd_louvain_multi_clustered_harmony_ws_41pcs_20240205_SCTonly_uvcm_comp_run2.rdata",
              "QCd_louvain_multi_clustered_harmony_ws_41pcs_20240216_HET_uvcm_v_abeta.rdata")

file_name <- c("40pcs_20240112_SCTonly_mgsubset_24hr_v_7d_abeta_comp",
               "41pcs_20240205_SCTonly_uvcm_comp_run2",
               "41pcs_20240216_HET_uvcm_v_abeta")
seurat_clustering <- c("0.3", "0.25", "0.4")

out_csv_folder <- "~/kzlinlab/projects/subject-de/out/kevin/Writeup4/"
out_plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup4/"

for(kk in 1:length(file_vec)){
  file <- file_vec[kk]
  print(paste0("Working on file: ", file))
  
  load(paste0(file_folder, file_vec[kk]))
  print("Done loading")
  
  Y <- table(ss_data_norm$Sample_ID, ss_data_norm@meta.data[,paste0("SCT_snn_res.", seurat_clustering[kk])])
  colnames(Y) <- paste0("cluster_", 0:(ncol(Y)-1))
  rowname_vec <- rownames(Y)
  Y <- as.data.frame(apply(as.matrix.noquote(Y),2,as.numeric))
  Y <- cbind(rowname_vec, Y)
  colnames(Y)[1] <- "SampleID"
  
  write.csv(Y, 
            row.names = FALSE,
            file = paste0(out_csv_folder, file_name[kk], "_composition.csv"))
  
  plot1 <- Seurat::DimPlot(ss_data_norm,
                           group.by = paste0("SCT_snn_res.", seurat_clustering[kk]))
  ggplot2::ggsave(filename = paste0(out_plot_folder, file_name[kk], "_umap.png"),
                  plot1, device = "png", width = 8, height = 5, units = "in",
                  dpi = 300)
}

print("Done! :)")