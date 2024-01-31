rm(list=ls())
library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
set.seed(10)

# remove the reference cells
keep_vec <- rep(FALSE, ncol(seurat_obj))
keep_vec[seurat_obj$CognitiveStatus != "Reference"] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# typical preprocessing
Seurat::DefaultAssay(seurat_obj) <- "RNA"
set.seed(10)
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj,
                             features = Seurat::VariableFeatures(object = seurat_obj),
                             verbose = F)
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              dims = 1:30)

#################

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
tab_mat <- table(seurat_obj$donor_id, seurat_obj$CognitiveStatus)
tab_mat <- tab_mat[,c("Dementia", "Nodementia")]

# now plot individuals
case_indiv <- rownames(tab_mat)[which(tab_mat[,"Dementia"] > 0)]
num_cases <- length(case_indiv)
case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(num_cases)
names(case_color_palette) <- case_indiv
control_indiv <- rownames(tab_mat)[which(tab_mat[,"Nodementia"] > 0)]
num_controls <- length(control_indiv)
control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(num_controls)
names(control_color_palette) <- control_indiv

col_palette <- c(case_color_palette, control_color_palette)

################

plot1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                         group.by = "donor_id",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup2/Writeup2_microglia-pvm_typical_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in",
                dpi = 300)
