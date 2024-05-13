rm(list=ls())
library(Seurat)

load("~/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_bad_raw-umap.RData")

keep_vec <- rep(FALSE, length(SeuratObject::Cells(ss_data_norm)))
names(keep_vec) <- SeuratObject::Cells(ss_data_norm)
keep_vec[rownames(metadata_df)] <- TRUE
ss_data_norm$keep <- keep_vec
ss_data_norm <- subset(ss_data_norm, keep == TRUE)
  
seurat_umap@cell.embeddings <- seurat_umap@cell.embeddings[SeuratObject::Cells(ss_data_norm),]
ss_data_norm[["umap"]] <- seurat_umap

var_vec <- c("Pt_ID", "Study_Designation", "SeqBatch")

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(ss_data_norm, 
                           reduction = "umap",
                           label = TRUE,
                           repel = TRUE,
                           label.size = 2.5,
                           group.by = variable)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("UMAP, Bad design, before integration"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  ggplot2::ggsave(plot1, 
                  file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_nature_bad_umap-", variable, ".png"),
                  height = 5, width = 8,
                  units = "in")
}




