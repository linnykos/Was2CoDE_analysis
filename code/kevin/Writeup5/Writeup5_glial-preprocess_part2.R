rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")

table(seurat_all$cell_type)
table(seurat_all$Subclass)
table(seurat_all$Supertype)
table(seurat_all$Cognitive.status)
table(seurat_all$Braak.stage)
table(seurat_all$Thal.phase)
table(seurat_all$CERAD.score)
table(seurat_all$APOE4.status)
table(seurat_all$LATE.NC.stage)
table(seurat_all$donor_id)
table(seurat_all$assay)
table(seurat_all$disease)
table(seurat_all$sex)

quantile(seurat_all$Continuous.Pseudo.progression.Score)

var_vec <- c("cell_type", "Subclass", "Supertype", "Cognitive.status",
             "Braak.stage", "Thal.phase", "CERAD.score", "APOE4.status",
             "LATE.NC.stage", "donor_id", "assay", "disease", "sex")

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_glial_umap-covariates.pdf",
    onefile = T, width = 8, height = 5)

for(variable in var_vec){
  plot1 <- Seurat::DimPlot(seurat_all, 
                           reduction = "umap",
                           group.by = variable,
                           raster = TRUE)
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  print(plot1)
}

plot1 <- Seurat::FeaturePlot(seurat_all,
                             features = "Continuous.Pseudo.progression.Score")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
print(plot1)

dev.off()
