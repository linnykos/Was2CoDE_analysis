rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_seurat-integration.RData")

var_vec <- c("Pt_ID", "Study_Designation", "SeqBatch")

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(ss_data_norm, 
                           reduction = "umap",
                           label = TRUE,
                           repel = TRUE,
                           label.size = 2.5,
                           group.by = variable)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("UMAP, Good design, after Seurat Integration"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  ggplot2::ggsave(plot1, 
                  file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_nature_good_seurat-integration_umap-", variable, ".png"),
                  height = 5, width = 8,
                  units = "in")
}
