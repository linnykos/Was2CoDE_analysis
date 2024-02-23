rm(list=ls())
library(Seurat)

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
count_mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, layer = "data", assay = "RNA"))
gene_vec <- intersect(colnames(esvd_mat), colnames(scvi_mat))

goi <- c("APOE", "PSEN1", "PSEN2", 
         "MYH11", "FZD3", "SORCS3",
         "TREM2", "CD33", "SPI1", "PU.1",
         "CR1")
goi[which(goi %in% gene_vec)]

mat_list <- list(
  count = count_mat[rownames(esvd_mat), gene_vec],
  esvd = esvd_mat[rownames(esvd_mat), gene_vec],
  scvi = scvi_mat[rownames(esvd_mat), gene_vec]
)

age_vec <- sapply(seurat_obj$development_stage, function(x){
  substr(x, start = 0, stop = 2)
})
age_vec[which(seurat_obj$development_stage == "80yearoldandoverhumanstage")] <- "90"
age_vec <- as.numeric(age_vec)
seurat_obj$AgeAtDeath <- age_vec

donor_vec <- droplevels(seurat_obj$donor_id[rownames(esvd_mat)])

violin_plot_func <- function(vec,
                             covariate_vec,
                             donor_vec,
                             bool_numeric = F){
  tab_mat <- table(donor_vec, covariate_vec)
  donors <- rownames(tab_mat)
  donors_covariate <- apply(tab_mat, 1, function(x){
    colnames(tab_mat)[which(x != 0)]
  })
  if(bool_numeric) donors_covariate <- as.numeric(donors_covariate)
  donors <- donors[order(donors_covariate)]
  donors_covariate <- donors_covariate[order(donors_covariate)]
  donors_alt <- sapply(1:length(donors_covariate), function(i){
    paste0(donors_covariate[i], ":", donors[i])
  })
  names(donors_alt) <- donors
  
  df <- data.frame(expression = vec,
                   donor = donors_alt[donor_vec])
  df_tmp <- df; df_tmp$donor <- as.factor(df_tmp$donor)
  anova_res <- stats::oneway.test(expression ~ donor, data = df_tmp)
  
  col_palette <- scales::hue_pal()(length(unique(donors_covariate)))
  names(col_palette) <- unique(donors_covariate)
  col_vec <- col_palette[as.character(donors_covariate)]
  names(col_vec) <- donors_alt
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(y=expression, x=donor))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=donor))
  p1 <- p1 + ggplot2::scale_fill_manual(values = col_vec) 
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::scale_x_discrete(limits = donors_alt,
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
    mat_list,
    bool_numeric = F
){
  vec_list <- lapply(mat_list, function(mat){
    vec <- mat[names(donor_vec), gene]
    scale(pmin(vec, stats::quantile(vec, probs = 0.99)))
  })
  
  p_list <- lapply(vec_list, function(vec){
    violin_plot_func(vec = vec,
                     covariate_vec = covariate_vec,
                     donor_vec = donor_vec,
                     bool_numeric = bool_numeric)
  })
  names(p_list) <- names(mat_list)
  
  p_list
}

####################

age_vec <- as.character(seurat_obj$AgeAtDeath[rownames(esvd_mat)])
gene <- "APOE"

p_list <- plot_multiple_violins(
  covariate_vec = age_vec,
  donor_vec = donor_vec,
  gene = gene,
  mat_list = mat_list,
  bool_numeric = T
)

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-donor_",
                                  gene, "_count.png"),
                p_list[["count"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-donor_",
                                  gene, "_eSVD.png"),
                p_list[["esvd"]], device = "png", width = 10, height = 3, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_violin-by-donor_",
                                  gene, "_scVI.png"),
                p_list[["scvi"]], device = "png", width = 10, height = 3, units = "in")

#############

sex_vec <- as.character(seurat_obj$sex[rownames(esvd_mat)])
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
                                  gene, "_scVI.png"),
                p_list[["scvi"]], device = "png", width = 10, height = 3, units = "in")

#############

sex_vec <- as.character(seurat_obj$sex[rownames(esvd_mat)])
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

education_vec <- as.character(seurat_obj$YearsOfEducation[rownames(esvd_mat)])
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
