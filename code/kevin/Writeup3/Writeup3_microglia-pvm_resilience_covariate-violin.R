rm(list=ls())
library(Seurat)


load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_preprocess_SCTransform.RData")
sct_mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                             layer = "data", 
                                             assay = "SCT"))

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
set.seed(10)

# scvi_mat <- read.csv("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_microglia-pvm_resilience_scvi_normalized-mat.csv")
# rownames(scvi_mat) <- scvi_mat[,1]
# scvi_mat <- scvi_mat[,-1]
# scvi_mat <- round(scvi_mat, 2)
# date_of_run <- Sys.time()
# session_info <- devtools::session_info()
# note <- "This is the scvi_normalized layer from the scVI fit of /home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia.h5ad"
# save(scvi_mat, date_of_run, session_info, note, file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_microglia-pvm_resilience_scvi_normalized-mat.RData")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_microglia-pvm_resilience_scvi_normalized-mat.RData")

esvd_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat) + tcrossprod(eSVD_obj$covariates[,"resilient_Resilient"], eSVD_obj$fit_Second$z_mat[,"resilient_Resilient"])
count_mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                               layer = "data", 
                                               assay = "RNA"))

gene_vec <- intersect(
  intersect(colnames(esvd_mat), colnames(scvi_mat)),
  colnames(sct_mat)
)

keep_vec <- rep(FALSE, length(SeuratObject::Cells(seurat_obj)))
keep_vec[which(SeuratObject::Cells(seurat_obj) %in% rownames(esvd_mat))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)
cell_vec <- SeuratObject::Cells(seurat_obj)

# https://alz-journals.onlinelibrary.wiley.com/doi/full/10.1002/alz.12181
goi <- c("APOE", "PSEN1", "PSEN2", 
         "MYH11", "FZD3", "SORCS3",
         "TREM2", "CD33", "SPI1", "PU.1",
         "CR1")
goi[which(goi %in% gene_vec)]

mat_list <- list(
  count = count_mat[cell_vec, gene_vec],
  esvd = esvd_mat[cell_vec, gene_vec],
  sct = sct_mat[cell_vec, gene_vec],
  scvi = as.matrix(log(scvi_mat[cell_vec, gene_vec]))
)

age_vec <- sapply(seurat_obj$development_stage, function(x){
  substr(x, start = 0, stop = 2)
})
age_vec[which(seurat_obj$development_stage == "80yearoldandoverhumanstage")] <- "90"
age_vec <- as.numeric(age_vec)
seurat_obj$AgeAtDeath <- age_vec

seurat_obj$donor_id <- droplevels(seurat_obj$donor_id)
donor_vec <- seurat_obj$donor_id

violin_plot_func <- function(vec,
                             covariate_vec,
                             donor_vec){
  tab_mat <- table(donor_vec, covariate_vec)
  donors <- rownames(tab_mat)
  donors_covariate <- apply(tab_mat, 1, function(x){
    colnames(tab_mat)[which(x != 0)]
  })
  names(donors_covariate) <- donors
  
  df <- data.frame(expression = vec,
                   donor = donor_vec)
  df_tmp <- df; df_tmp$donor <- as.factor(df_tmp$donor)
  anova_res <- stats::oneway.test(expression ~ donor, data = df_tmp)
  
  col_palette <- scales::hue_pal()(length(unique(covariate_vec)))
  names(col_palette) <- unique(covariate_vec)
  col_vec <- col_palette[donors_covariate]
  names(col_vec) <- donors
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(y=expression, x=donor))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=donor))
  p1 <- p1 + ggplot2::scale_fill_manual(values = col_vec) 
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::scale_x_discrete(limits = donors[order(donors_covariate)],
                                       guide = ggplot2::guide_axis(angle = 45))
  p1 <- p1 + ggplot2::ylab("Expression")
  p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  p1 <- p1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2)))
  p1
}

plot_multiple_violins <- function(
    covariate_vec,
    donor_vec,
    gene,
    mat_list
){
  vec_list <- lapply(mat_list, function(mat){
    vec <- mat[, gene]
    pmax(vec, pmin(vec, stats::quantile(vec, probs = 0.95)), stats::quantile(vec, probs = 0.05))
  })
  
  p_list <- lapply(vec_list, function(vec){
    violin_plot_func(vec = vec,
                     covariate_vec = covariate_vec,
                     donor_vec = donor_vec)
  })
  names(p_list) <- names(mat_list)
  
  p_list
}

####################

age_vec <- as.character(seurat_obj$AgeAtDeath)
gene <- "APOE"

p_list <- plot_multiple_violins(
  covariate_vec = age_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-donor_",
                                  gene, "_count.png"),
                p_list[["count"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-donor_",
                                  gene, "_eSVD.png"),
                p_list[["esvd"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-donor_",
                                  gene, "_SCTransform.png"),
                p_list[["sct"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-donor_",
                                  gene, "_scVI.png"),
                p_list[["scvi"]], device = "png", width = 10, height = 3, units = "in")

#############

sex_vec <- as.character(seurat_obj$sex)
gene <- "APOE"

p_list <- plot_multiple_violins(
  covariate_vec = sex_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_count.png"),
                p_list[["count"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_eSVD.png"),
                p_list[["esvd"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_SCTransform.png"),
                p_list[["sct"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_scVI.png"),
                p_list[["scvi"]], device = "png", width = 10, height = 3, units = "in")

#############

sex_vec <- as.character(seurat_obj$sex)
gene <- "CR1" # "APOE"

p_list <- plot_multiple_violins(
  covariate_vec = sex_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_count.png"),
                p_list[["count"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_eSVD.png"),
                p_list[["esvd"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_scVI.png"),
                p_list[["scvi"]], device = "png", width = 10, height = 3, units = "in")

############

education_vec <- as.character(seurat_obj$YearsOfEducation)
gene <- "APOE"

p_list <- plot_multiple_violins(
  covariate_vec = education_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-education_",
                                  gene, "_count.png"),
                p_list[["count"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-education_",
                                  gene, "_eSVD.png"),
                p_list[["esvd"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-education_",
                                  gene, "_scVI.png"),
                p_list[["scvi"]], device = "png", width = 10, height = 3, units = "in")

#####################################

sex_vec <- seurat_obj$sex
donor_vec <- seurat_obj$donor_id

tab_mat <- table(donor_vec, sex_vec)
donors <- rownames(tab_mat)
donors_covariate <- apply(tab_mat, 1, function(x){
  colnames(tab_mat)[which(x != 0)]
})
names(donors_covariate) <- donors

donor_list <- lapply(donors, function(donor){
  which(donor_vec == donor)
})
median_expression_list <- lapply(mat_list, function(mat){
  tmp <- sapply(donor_list, function(idx_vec){
    if(is.matrix(mat)){
      matrixStats::colMedians(mat[idx_vec,])
    } else {
      sparseMatrixStats::colMedians(mat[idx_vec,])
    }
  })
  rownames(tmp) <- gene_vec
  colnames(tmp) <- donors
  tmp
})
names(median_expression_list) <- names(mat_list)

esvd_vec <- sapply(1:p, function(j){
  tmp <- stats::t.test(
    x = median_expression_list[["esvd"]][j,donors_covariate == "female"],
    y = median_expression_list[["esvd"]][j,donors_covariate == "male"]
  )
  -log10(tmp$p.value)
})
scvi_vec <- sapply(1:p, function(j){
  vec1 <- mat_list[["scvi"]][f_idx,j]
  vec2 <- mat_list[["scvi"]][m_idx,j]
  vec1[is.infinite(vec1)] <- -5
  vec2[is.infinite(vec2)] <- -5
  
  tmp <- stats::t.test(
    x = vec1,
    y = vec2,
  )
  -log10(tmp$p.value)
})
sct_vec <- sapply(1:p, function(j){
  tmp <- stats::t.test(
    x = mat_list[["sct"]][f_idx,j],
    y = mat_list[["sct"]][m_idx,j],
  )
  -log10(tmp$p.value)
})
sct_vec[is.na(sct_vec)] <- 0


round(quantile(sct_vec),2); 
round(quantile(scvi_vec),2); 
round(quantile(esvd_vec),2)

idx <- intersect(
  intersect(
    which(sct_vec > 4),
    which(scvi_vec > 50)
  ),
  which(esvd_vec < 0.2)
)
length(idx)
gene_vec[idx]

###########

sex_vec <- as.character(seurat_obj$sex)
gene <- "TLE4" 

p_list <- plot_multiple_violins(
  covariate_vec = sex_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_eSVD.png"),
                p_list[["esvd"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_scVI.png"),
                p_list[["scvi"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_SCTransform.png"),
                p_list[["sct"]], device = "png", width = 10, height = 3, units = "in")


