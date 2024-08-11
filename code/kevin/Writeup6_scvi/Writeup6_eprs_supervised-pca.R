rm(list=ls())
library(Seurat)
source("~/kzlinlab/projects/subject-de/git/subject-de_kevin/code/kevin/Writeup6_scvi/supervised_pca.R")

load("~/kzlinlab/data/jayadev_pu1-only_eprs/QCd_norm_anchored_rpca__ws_anch_mg_subset.rdata")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
ss_data_norm <- Seurat::UpdateSeuratObject(ss_data_norm)

colnames(ss_data_norm@meta.data)
table(ss_data_norm$UWA, ss_data_norm$ePRS)

scvi_mat <- read.csv("~/kzlinlab/projects/subject-de/git/subject-de_kevin/csv/kevin/Writeup6/Writeup6_eprs_mg_scvi_scVI.csv")

rownames(scvi_mat) <- scvi_mat[,1]
scvi_mat <- scvi_mat[,-1]
scvi_mat <- as.matrix(scvi_mat)
scvi_mat <- scvi_mat[Seurat::Cells(ss_data_norm),]
colnames(scvi_mat) <- paste0("scVI_", 1:ncol(scvi_mat))

ss_data_norm[["scVI"]] <- Seurat::CreateDimReducObject(embeddings = scvi_mat, 
                                                       assay = "integrated")

#########

# now for supervised PCA (for 1D)

y <- matrix(1, nrow = length(Seurat::Cells(ss_data_norm)), ncol = 3)
metadata <- ss_data_norm@meta.data

y[which(metadata$ePRS == "High"),1] <- -1
y[which(metadata$sex == "Male"),2] <- -1
y[which(metadata$cognitive_status == "dementia"),3] <- -1

spca_res <- supervised_pca(
  x = scvi_mat,
  y = y,
  k = 1
)

ss_data_norm[["sPCA"]] <- Seurat::CreateDimReducObject(embeddings = spca_res$dimred, 
                                                       assay = "integrated")

plot1 <- Seurat::DimPlot(ss_data_norm, 
                         reduction = "sPCA",
                         group.by = "ePRS",
                         raster = TRUE)
ggplot2::ggsave(plot1,
                filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_eprs_mg_scvi_sPCA_ePRS.png",
                width = 8, height = 8)

############

# now for supervised PCA (visualization)

y <- matrix(1, nrow = length(Seurat::Cells(ss_data_norm)), ncol = 3)
metadata <- ss_data_norm@meta.data

y[which(metadata$ePRS == "High"),1] <- -1
y[which(metadata$sex == "Male"),2] <- -1
y[which(metadata$cognitive_status == "dementia"),3] <- -1

spca_res <- supervised_pca(
  x = scvi_mat,
  y = y
)

ss_data_norm[["sPCA"]] <- Seurat::CreateDimReducObject(embeddings = spca_res$dimred, 
                                                       assay = "integrated")

plot1 <- Seurat::DimPlot(ss_data_norm, 
                         reduction = "sPCA",
                         group.by = "ePRS",
                         raster = TRUE)
ggplot2::ggsave(plot1,
                filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_eprs_mg_scvi_sPCA_ePRS.png",
                width = 8, height = 8)

##################

# try a different strategy

spca_res <- supervised_pca(
  x = scvi_mat,
  y = y,
  k = 10
)

set.seed(10)
umap_res <- Seurat::RunUMAP(spca_res$dimred)
ss_data_norm[["sPCA2"]] <- umap_res

plot1 <- Seurat::DimPlot(ss_data_norm, 
                         reduction = "sPCA2",
                         group.by = "ePRS",
                         raster = TRUE)
ggplot2::ggsave(plot1,
                filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_eprs_mg_scvi_sPCA-v2_ePRS.png",
                width = 8, height = 8)


###############

# turns out we don't need supervised-pca

set.seed(10)
umap_res <- Seurat::RunUMAP(scvi_mat)
ss_data_norm[["scVI_UMAP"]] <- umap_res

plot1 <- Seurat::DimPlot(ss_data_norm, 
                         reduction = "scVI_UMAP",
                         group.by = "ePRS",
                         raster = TRUE)
ggplot2::ggsave(plot1,
                filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_eprs_mg_scvi_ePRS_R-umap.png",
                width = 8, height = 8)

umap_res <- ss_data_norm[["scVI_UMAP"]]@cell.embeddings
n <- nrow(umap_res)
color_vec <- rep("darkgreen", n)
idx <- which(ss_data_norm$ePRS == "High")
color_vec[idx] <- "orange"

set.seed(10)
order_vec <- sample(n)

png("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_eprs_mg_scvi_ePRS_R-umap_full.png",
    height = 1200, width = 1500, units = "px", res = 300)
par(mar = c(2, 2, 0.5, 0.5))
plot(umap_res[order_vec,1],
     umap_res[order_vec,2],
     pch = 16,
     cex = 0.2,
     col = color_vec[order_vec],
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     bty = "n")
axis(1); axis(2)
graphics.off()


