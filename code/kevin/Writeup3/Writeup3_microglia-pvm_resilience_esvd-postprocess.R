rm(list=ls())
library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
set.seed(10)

# remove the reference cells
keep_vec <- rep(FALSE, length(SeuratObject::Cells(seurat_obj)))
names(keep_vec) <- SeuratObject::Cells(seurat_obj)
keep_vec[rownames(eSVD_obj$dat)] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

x <- scale(eSVD_obj$fit_Second$x_mat)
set.seed(10)
umap_esvd <- Seurat::RunUMAP(x)
seurat_obj[["esvd"]] <- Seurat::CreateDimReducObject(
  embeddings = umap_esvd@cell.embeddings[SeuratObject::Cells(seurat_obj),]
)

##########################

donor_vec <- seurat_obj$donor_id
names(donor_vec) <- SeuratObject::Cells(seurat_obj)
donor_vec <- donor_vec[rownames(eSVD_obj$dat)]
donor_vec <- droplevels(donor_vec)
resilient_vec <- eSVD_obj$covariates[,"resilient_Resilient"]

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
tab_mat <- table(donor_vec, resilient_vec)
tab_mat <- tab_mat[,c("0", "1")]

# now plot individuals
case_indiv <- rownames(tab_mat)[which(tab_mat[,"0"] > 0)]
num_cases <- length(case_indiv)
case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(num_cases)
names(case_color_palette) <- case_indiv
control_indiv <- rownames(tab_mat)[which(tab_mat[,"1"] > 0)]
num_controls <- length(control_indiv)
control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(num_controls)
names(control_color_palette) <- control_indiv

col_palette <- c(case_color_palette, control_color_palette)

plot1 <- Seurat::DimPlot(seurat_obj, reduction = "esvd",
                         group.by = "donor_id",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_microglia-pvm_esvd_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in",
                dpi = 300)
