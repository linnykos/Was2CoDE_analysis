rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_glial_tsne.RData")

source("Writeup5_glial_palette.R")

supertype_vec <- seurat_all$Supertype
supertype_vec[supertype_vec == "Astro_6-SEAAD"] <- "Astro_6"
supertype_vec[supertype_vec == "Micro-PVM_1_1-SEAAD"] <- "Micro-PVM_1"
supertype_vec[supertype_vec == "Micro-PVM_2_1-SEAAD"] <- "Micro-PVM_2"
supertype_vec[supertype_vec == "Micro-PVM_2_2-SEAAD"] <- "Micro-PVM_2"
supertype_vec[supertype_vec == "Micro-PVM_2_3-SEAAD"] <- "Micro-PVM_2"
supertype_vec[supertype_vec == "Micro-PVM_3-SEAAD"] <- "Micro-PVM_3"
supertype_vec[supertype_vec == "Micro-PVM_4-SEAAD"] <- "Micro-PVM_4"
supertype_vec[supertype_vec == "Oligo_2_1-SEAAD"] <- "Oligo_2"
supertype_vec[supertype_vec == "OPC_2_1-SEAAD"] <- "OPC_2"
supertype_vec[supertype_vec == "OPC_2_2-SEAAD"] <- "OPC_2"
seurat_all$Supertype2 <- supertype_vec

table(seurat_all$Supertype2)

for(i in 1:length(tsne_list)){
  print(paste("Working on trial:", i))
  set.seed(i)
  seurat_all[["tsne"]] <- tsne_list[[i]]
  
  plot1 <- Seurat::DimPlot(seurat_all, 
                           reduction = "tsne",
                           label = TRUE,
                           repel = TRUE,
                           label.size = 2.5,
                           group.by = "Supertype2",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("t-SNE, Seed: ", i))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  ggplot2::ggsave(plot1, 
                  file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_glial_tsne-", i, ".png"),
                  height = 5, width = 5,
                  units = "in")
}
