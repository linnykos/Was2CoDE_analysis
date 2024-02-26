rm(list=ls())
library(Seurat)


load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_preprocess_SCTransform.RData")
sct_mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                             layer = "data", 
                                             assay = "SCT"))

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
set.seed(10)


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
seurat_obj2 <- seurat_obj

# subset the cells
donor_ids <- c("H2033033", "H2133005", "H2133020", "H2133046",
               "H2133007", "H2033018", "H2033045", "H2133036")
keep_vec <- rep(FALSE, length(SeuratObject::Cells(seurat_obj)))
keep_vec[which(seurat_obj$donor_id %in% donor_ids)] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)
cell_vec <- SeuratObject::Cells(seurat_obj)
seurat_obj$donor_id <- droplevels(seurat_obj$donor_id)
donor_vec <- seurat_obj$donor_id

mat_list <- list(
  esvd = esvd_mat[cell_vec, gene_vec],
  sct = sct_mat[cell_vec, gene_vec],
  scvi = as.matrix(log(scvi_mat[cell_vec, gene_vec]))
)

violin_plot_func <- function(vec,
                             covariate_vec,
                             donor_vec){
  vec <- vec - min(vec)
  
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
  
  quant_vec <- stats::quantile(vec, probs = c(0.05, 0.95))
  p1 <- ggplot2::ggplot(df, ggplot2::aes(y=expression, x=donor))
  p1 <- p1 + ggplot2::geom_boxplot(ggplot2::aes(fill=donor))
  p1 <- p1 + ggplot2::scale_fill_manual(values = col_vec) 
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::ylim(quant_vec[1], quant_vec[2])
  p1 <- p1 + ggplot2::scale_x_discrete(limits = donors[order(donors_covariate)],
                                       guide = ggplot2::guide_axis(angle = 45))
  p1 <- p1 + ggplot2::ylab("Expression")
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
  })
  
  p_list <- lapply(1:length(vec_list), function(ii){
    violin_plot_func(vec = vec_list[[ii]],
                     covariate_vec = covariate_vec,
                     donor_vec = donor_vec)
  })
  names(p_list) <- names(mat_list)
  
  p_list
}

#################

sex_vec <- as.character(seurat_obj$sex)
gene <- "APOE"

p_list <- plot_multiple_violins(
  covariate_vec = sex_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_eSVD_short.png"),
                p_list[["esvd"]], device = "png", width = 4, height = 2, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_SCTransform_short.png"),
                p_list[["sct"]], device = "png", width = 4, height = 2, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_scVI_short.png"),
                p_list[["scvi"]], device = "png", width = 4, height = 2, units = "in")

#################

sex_vec <- as.character(seurat_obj$sex)
gene <- "TTTY14"

p_list <- plot_multiple_violins(
  covariate_vec = sex_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_eSVD_short.png"),
                p_list[["esvd"]], device = "png", width = 4, height = 2, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_SCTransform_short.png"),
                p_list[["sct"]], device = "png", width = 4, height = 2, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-sex_",
                                  gene, "_scVI_short.png"),
                p_list[["scvi"]], device = "png", width = 4, height = 2, units = "in")

###########

scvi_vec <- mat_list[["scvi"]][,"APOE"]
names(scvi_vec) <- rownames(mat_list[["scvi"]])
donor_vec <- seurat_obj$donor_id
mean(scvi_vec[which(donor_vec=="H2033008")])
