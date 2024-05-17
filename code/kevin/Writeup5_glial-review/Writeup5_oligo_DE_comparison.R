rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")

seurat_all <- subset(seurat_all, Supertype == "Oligo_2")

quantile(seurat_all$Continuous.Pseudo.progression.Score)
table(seurat_all$Continuous.Pseudo.progression.Score, seurat_all$donor_id)

median_val <- stats::median(seurat_all$Continuous.Pseudo.progression.Score)
donor_vec <- rep("low", length(Seurat::Cells(seurat_all)))
donor_vec[seurat_all$Continuous.Pseudo.progression.Score > median_val] <- "high"
seurat_all$classification <- as.factor(donor_vec)

set.seed(10)
Seurat::Idents(seurat_all) <- "classification"
de_markers_good <- Seurat::FindMarkers(seurat_all, 
                                       ident.1 = "low", 
                                       ident.2 = "high",
                                       min.pct = 0,
                                       logfc.threshold = 0)

##########

pca_obj <- seurat_all[["pca"]]
num_pc <- 30
recon_data <- tcrossprod(pca_obj@cell.embeddings[,1:num_pc], pca_obj@feature.loadings[,1:num_pc])

# Seurat::FindMarkers doesn't like negative numbers
recon_data <- recon_data - min(recon_data) + 1

SeuratObject::LayerData(seurat_all, 
                        layer = "data",
                        assay = "RNA") <- t(recon_data)

set.seed(10)
Seurat::Idents(seurat_all) <- "classification"
de_markers_bad <- Seurat::FindMarkers(seurat_all, 
                                  assay = "RNA",
                                  layer = "counts",
                                  ident.1 = "low", 
                                  ident.2 = "high",
                                  features = SeuratObject::Features(seurat_all),
                                  min.pct = 0,
                                  logfc.threshold = 0)

de_markers_bad <- de_markers_bad[rownames(de_markers_good),]

tmp <- pmax(de_markers_bad$p_val, 10e-300)
df <- data.frame(good = -log10(de_markers_good$p_val),
                 bad = -log10(tmp))

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=bad, y=good)) 
plot1 <- plot1 + ggplot2::geom_point()
plot1 <- plot1 + ggplot2::labs(x = "Improper DE (-Log10 p-value)",
                               y = "Proper DE (-Log10 p-value)",
                               tilte = paste0("Correlation: ", round(stats::cor(df[,1], df[,2]), 2)))
plot1 <- plot1 + ggplot2::coord_fixed()

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/Writeup5_oligo_comparison.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")



