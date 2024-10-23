rm(list=ls())
library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_preprocess_SCTransform.RData")
sct_mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                             layer = "data", 
                                             assay = "SCT"))

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
esvd_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat) + tcrossprod(eSVD_obj$covariates[,"resilient_Resilient"], eSVD_obj$fit_Second$z_mat[,"resilient_Resilient"])

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_microglia-pvm_resilience_scvi_normalized-mat.RData")
count_mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                               layer = "data", 
                                               assay = "RNA"))

set.seed(10)

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/"

###############
scvi_mat <- as.matrix(scvi_mat)

gene_vec <- intersect(colnames(esvd_mat), colnames(scvi_mat))
gene_vec <- intersect(gene_vec, colnames(sct_mat))
cell_names <- rownames(sct_mat)

sct_mat <- sct_mat[cell_names,gene_vec]
esvd_mat <- esvd_mat[cell_names,gene_vec]
scvi_mat <- scvi_mat[cell_names,gene_vec]
count_mat <- count_mat[cell_names,gene_vec]

# compute a bunch of PCAs
set.seed(10); sct_pca <- irlba::irlba(sct_mat, nv = 10)$u
set.seed(10); esvd_pca <- irlba::irlba(esvd_mat, nv = 10)$u
set.seed(10); scvi_pca <- irlba::irlba(scvi_mat, nv = 10)$u
set.seed(10); count_pca <- irlba::irlba(count_mat, nv = 10)$u

# get the sex vector
sex_vec <- seurat_obj@meta.data[cell_names,"sex"]
sex_vec <- as.numeric(factor(sex_vec))-1

pca_list <- list(sct = sct_pca,
                 esvd = esvd_pca,
                 scvi = scvi_pca,
                 count = count_pca)

for(kk in 1:length(pca_list)){
  print(names(pca_list)[kk])
  pca_mat <- pca_list[[kk]]
  data <- data.frame(cbind(sex_vec, pca_mat))
  colnames(data)[1] <- "sex"
  logit_res <- stats::glm(sex ~ ., data = data, family = "binomial")
  predict_res <- stats::predict(logit_res, newdata=data, type="response")
  print(table(sex_vec, predict_res <= 0.5))
  print("===")
}

# compute the ROC curves
roc_func <- function(true_vec, prob_vec){
  stopifnot(all(prob_vec >= 0), all(prob_vec <= 1),
            all(true_vec %in% c(0,1)))
  seq_val <- seq(0, 1, length.out = 100)
  true_1_idx <- which(true_vec == 1)
  true_0_idx <- which(true_vec == 0)
  
  tpr_vec <- sapply(seq_val, function(x){
    bool_vec <- (prob_vec <= x)
    length(intersect(true_1_idx, which(bool_vec == 1)))/length(true_1_idx)
  })
  fpr_vec <- sapply(seq_val, function(x){
    bool_vec <- (prob_vec <= x)
    length(intersect(true_0_idx, which(bool_vec == 1)))/length(true_0_idx)
  })
  
  cbind(seq_val, tpr_vec, fpr_vec)
}

roc_list <- lapply(1:length(pca_list), function(kk){
  pca_mat <- pca_list[[kk]]
  data <- data.frame(cbind(sex_vec, pca_mat))
  colnames(data)[1] <- "sex"
  logit_res <- stats::glm(sex ~ ., data = data, family = "binomial")
  predict_res <- stats::predict(logit_res, newdata=data, type="response")
  roc_func(true_vec = sex_vec, 
           prob_vec = predict_res)
})
names(roc_list) <- names(pca_list)

col_vec <- c(rgb(122, 49, 126, maxColorValue = 255),
             rgb(236, 134, 46, maxColorValue = 255),
             rgb(50, 174, 255, maxColorValue = 255),
             "black")
names(col_vec) <- names(pca_list)

png(paste0(plot_folder, "Writeup3_sex_ROC.png"),
    height = 1000, width = 800, res = 300, units = "px")
plot(NA, 
     xlim = c(0,1), 
     ylim = c(0,1),
     xlab = "False positive rate",
     ylab = "True positive rate",
     xaxt = "n",
     yaxt = "n",
     bty = "n",
     asp = TRUE)
axis(1); axis(2)
lines(x = c(0,1), y = c(0,1), lty = 2)
for(kk in 1:length(roc_list)){
  lines(y = roc_list[[kk]][,3],
         x = roc_list[[kk]][,2],
         lwd = 2.5,
         col = col_vec[kk])
}
dev.off()
