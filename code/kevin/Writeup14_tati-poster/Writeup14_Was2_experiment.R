rm(list=ls())

library(Seurat)

df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/out/tati/Writeup5/Prater_dataset_ingredients.csv")
rownames(df) <- df$Gene
load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/out/tati/Writeup5/Writeup5_microglia_ideascustom.RData")

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

# compute the logFC
dat_mat <- SeuratObject::LayerData(
  object = ss_data_norm,
  layer = "data",
  assay = "RNA"
)
case_idx <- which(ss_data_norm$Study_Designation == "AD")
control_idx <- which(ss_data_norm$Study_Designation == "Ctrl")

case_mean <- Matrix::rowMeans(dat_mat[,case_idx])
control_mean <- Matrix::rowMeans(dat_mat[,control_idx])
diff_vec <- case_mean - control_mean
names(diff_vec) <- rownames(dat_mat)

deseq2_lfc <- df[,"DESeq2_logFC"]
names(deseq2_lfc) <- df$Gene

diff_vec <- diff_vec[names(deseq2_lfc)]
diff_vec <- rank(diff_vec)
deseq2_lfc <- rank(deseq2_lfc)

cor_val <- stats::cor(deseq2_lfc, diff_vec)
plot(deseq2_lfc, diff_vec, main = paste0("Cor: ", round(cor_val, 2)), 
     pch = 16,
     col = rgb(0.5,0.5,0.5,0.5))
