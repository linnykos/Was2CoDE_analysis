rm(list=ls())
library(Seurat)

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup14/"
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

table(ss_data_norm$Study_Designation)
color_vec <- c(AD = "#6a0dad", Ctrl = "#f4c20d")

plot1 <- Seurat::DimPlot(ss_data_norm, 
                         reduction = "umap",
                         group.by = "Study_Designation",
                         cols = color_vec,
                         raster = TRUE,
                         raster.dpi = c(1000, 1000),
                         pt.size = 2)
plot1 <- plot1 + Seurat::NoLegend() + ggplot2::labs(x = "", y = "", title = "") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )

ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup14_prater-umap.pdf"),
                height = 4, width = 4)

#######

out_folder <- "~/kzlinlab/projects/subject-de/out/kevin/Writeup10/"
load(paste0(out_folder, "Writeup10_sea-ad_microglia_scVI-postprocessed.RData"))

table(seurat_obj$ADNC)
color_vec <- c(Case = "#6a0dad", Control = "#f4c20d")

plot1 <- Seurat::DimPlot(seurat_obj, 
                         reduction = "umap",
                         group.by = "ADNC",
                         cols = color_vec,
                         raster = TRUE,
                         raster.dpi = c(1000, 1000),
                         pt.size = 2)
plot1 <- plot1 + Seurat::NoLegend() + ggplot2::labs(x = "", y = "", title = "") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )

ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup14_seaad-umap.pdf"),
                height = 4, width = 4)


#######

out_folder <- "~/kzlinlab/projects/subject-de/out/kevin/Writeup11/"
load(paste0(out_folder, "Writeup11_rosmap_scVI-postprocessed.RData"))

table(seurat_obj$ADpath)
color_vec <- c(yes = "#6a0dad", no = "#f4c20d")

plot1 <- Seurat::DimPlot(seurat_obj, 
                         reduction = "umap",
                         group.by = "ADpath",
                         cols = color_vec,
                         raster = TRUE,
                         raster.dpi = c(1000, 1000),
                         pt.size = 2)
plot1 <- plot1 + Seurat::NoLegend() + ggplot2::labs(x = "", y = "", title = "") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  )

ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup14_rosmap-umap.pdf"),
                height = 4, width = 4)

################

df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Prater_dataset_ingredients.csv")

