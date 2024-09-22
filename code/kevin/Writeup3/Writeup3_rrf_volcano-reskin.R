rm(list=ls())
library(EnhancedVolcano)
library(ggplot2)
set.seed(10)

out_folder <- "~/kzlinlab/projects/fxomics_2024_kevin/out/kevin/Writeup2/"
plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/"
file_prefix <- "Writeup2_pu1-20human-highADNC-sorted_cluster-c:"
file_suffix <- "_eSVD.RData"
cluster_idx <- 10

load(paste0(out_folder, file_prefix, cluster_idx, file_suffix))

lfc_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
pvalue_vec <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
pval_adj_vec <- eSVD_obj$pvalue_list$fdr_vec
idx <- which(pval_adj_vec <= 0.05)
pCutoff <- max(pvalue_vec[idx])
FCcutoff <- 0
xlim <- c(-1,1) * quantile(abs(lfc_vec), probs = 0.99)

idx <- which(abs(pvalue_vec) <= max(xlim))
ylim <- c(0, max(-log10(pvalue_vec[idx])))

res <- data.frame(gene = names(lfc_vec),
                  avg_log2FC = lfc_vec,
                  p_val = pvalue_vec,
                  p_val_adj = pval_adj_vec)
res <- res[order(res$p_val, decreasing = FALSE),]

plot1 <- EnhancedVolcano::EnhancedVolcano(
  res,
  lab = res$gene,
  x = "avg_log2FC",
  y = "p_val",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  col = c('grey30', 'grey30', 'grey30', 'red2'),
  xlim = xlim,
  ylim = ylim
)

ggplot2::ggsave(plot1,
                filename = paste0(plot_folder, "Writeup3_rrf_volcano-reskin.png"),
                width = 5.5, height = 7)
