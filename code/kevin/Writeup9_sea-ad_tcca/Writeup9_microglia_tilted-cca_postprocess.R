# https://linnykos.github.io/tiltedCCA/articles/embryo.html
rm(list=ls())

library(Seurat)
library(tiltedCCA)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_tcca.RData")

# first make plots
categorical_vec <- c("Cognitive status",
                     "APOE4 status",
                     "LATE-NC stage", "CERAD score",
                     "RNA_snn_res.0.5", 
                     "ATAC_snn_res.0.5")
plot_filename <- c("tcca-common", "tcca-RNA-distinct", "tcca-ATAC-distinct")
umap_name <- c("common_tcca", "distinct1_tcca", "distinct2_tcca")

for(kk in 1:length(plot_filename)){
  print(plot_filename[kk])
  
  file <- plot_filename[kk]
  umap <- umap_name[kk]
  
  pdf(paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup9/Writeup9_microglia_", file,"_plots.pdf"),
      onefile = T, width = 7, height = 5)
  
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = umap,
                           group.by = "donor_id", 
                           cols = col_vec)
  print(plot1)
  
  for(variable in categorical_vec){
    plot1 <- Seurat::DimPlot(seurat_obj, 
                             reduction = umap,
                             group.by = variable)
    print(plot1)
  }
  
  dev.off()
}

# now for synchrony score
# from https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_steadystate.R
rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
seurat_obj$synchrony <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(seurat_obj$synchrony), max(seurat_obj$synchrony), length.out = num_color)
color_vec <- sapply(seurat_obj$synchrony, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

names(alignment_vec) <- Seurat::Cells(seurat_obj)


png(paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup9/Writeup9_microglia_synchrony.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = seurat_obj[["common_tcca"]]@cell.embeddings[,1],
     y = seurat_obj[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(seurat_obj[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(seurat_obj[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Alignment between ATAC and common RNA\nYellow-high, Green-low"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()

###########

# make a histogram of pseudoprogression by alignment

df <- data.frame(
  donor_id = seurat_obj$donor_id,
  pseudoprogression = seurat_obj@meta.data[,"Continuous Pseudo-progression Score"],
  synchrony = seurat_obj$synchrony,
  cognitive_status = seurat_obj@meta.data[,"Cognitive status"]
)

# find the ordering of donors based on pseudoprogression
donor_vec <- unique(seurat_obj$donor_id)
pseudoprogression_vec <- sapply(donor_vec, function(donor){
  idx <- which(df$donor_id == donor)
  mean(df$pseudoprogression[idx])
})
donor_vec <- donor_vec[order(pseudoprogression_vec, decreasing = FALSE)]

df$donor_id <- factor(df$donor_id, levels = donor_vec)
df$cognitive_status <- factor(df$cognitive_status)

col_vec <- c("Dementia" = "red", "No dementia" = "blue")
cor_val <- stats::cor(df$synchrony, df$pseudoprogression)

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = donor_id, y = synchrony))
p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill = cognitive_status))
p1 <- p1 + ggplot2::scale_fill_manual(values = col_vec) 
p1 <- p1 + ggplot2::geom_jitter(shape = 16, 
                                position = ggplot2::position_jitter(0.2), 
                                alpha = 0.3, 
                                size = 0.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
p1 <- p1 + ggplot2::scale_x_discrete(limits = levels(df$donor_id),
                                     guide = ggplot2::guide_axis(angle = 45))
p1 <- p1 + ggplot2::ylab("Synchrony score")
p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="white")
p1 <- p1 + ggplot2::ggtitle(paste("Correlation:", round(cor_val, 2)))
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup9/Writeup9_microglia_synchrony-violinplots.png"),
                p1, device = "png", width = 12, height = 3, units = "in")

