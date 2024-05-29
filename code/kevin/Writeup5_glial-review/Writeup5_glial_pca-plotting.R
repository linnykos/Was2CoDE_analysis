rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_glial_umaps.RData")

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

plot1 <- Seurat::DimPlot(seurat_all, 
                         reduction = "pca",
                         dims = c(1,2),
                         label = TRUE,
                         repel = TRUE,
                         label.size = 2.5,
                         group.by = "Supertype2",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("PCA, 1 vs. 2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

ggplot2::ggsave(plot1, 
                file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_glial_pca_1-2.png"),
                height = 5, width = 5,
                units = "in")

##########

plot1 <- Seurat::DimPlot(seurat_all, 
                         reduction = "pca",
                         dims = c(1,3),
                         label = TRUE,
                         repel = TRUE,
                         label.size = 2.5,
                         group.by = "Supertype2",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("PCA, 1 vs. 3"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

ggplot2::ggsave(plot1, 
                file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_glial_pca_1-3.png"),
                height = 5, width = 5,
                units = "in")

##########

plot1 <- Seurat::DimPlot(seurat_all, 
                         reduction = "pca",
                         dims = c(2,3),
                         label = TRUE,
                         repel = TRUE,
                         label.size = 2.5,
                         group.by = "Supertype2",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("PCA, 2 vs. 3"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

ggplot2::ggsave(plot1, 
                file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_glial_pca_2-3.png"),
                height = 5, width = 5,
                units = "in")

##########

plot1 <- Seurat::DimPlot(seurat_all, 
                         reduction = "pca",
                         dims = c(3,4),
                         label = TRUE,
                         repel = TRUE,
                         label.size = 2.5,
                         group.by = "Supertype2",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("PCA, 3 vs. 4"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

ggplot2::ggsave(plot1, 
                file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_glial_pca_3-4.png"),
                height = 5, width = 5,
                units = "in")