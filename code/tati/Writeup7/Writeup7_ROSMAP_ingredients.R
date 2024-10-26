rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(IdeasCustom)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(dplyr)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_NEBULA.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_Pseudobulk-DEseq2.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_eSVD.RData")
# load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_was2_wilcox.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_variance_results.RData")


dat_list <- list(
  eSVD = list(
    pvalue = 10^(-eSVD_obj$pvalue_list$log10pvalue),
    logFC = log2(eSVD_obj$case_mean / eSVD_obj$control),
    genes = names(eSVD_obj$teststat_vec)
  ),
  # WAS2 = list(
  #   pvalue = results_mat[,"p_val"],
  #   logFC = log2((results_mat[, "mean_dn"]) / ((results_mat[, "mean_nn"] + results_mat[, "mean_dd"]) / 2)),
  #   genes = rownames(results_mat)
  # ),
  DESeq2 = list(
    pvalue = deseq2_res[,"pvalue"],
    logFC = deseq2_res[,"log2FoldChange"],
    genes = rownames(deseq2_res)
  ),
  NEBULA = list(
    pvalue = nebula_res$summary[,"p_ADpathyes"],
    logFC = nebula_res$summary[,"logFC_ADpathyes"],
    genes = nebula_res$summary$gene
  )
)

all_genes <- unique(c(
  dat_list$eSVD$genes,
  rownames(deseq2_res),
  dat_list$NEBULA$genes,
  # dat_list$WAS2$genes,
  result_df$Gene
))

combined_df <- data.frame(Gene = all_genes)

#safely match and fill data
safe_match <- function(data, genes, all_genes) {
  matched <- match(all_genes, genes)
  ifelse(is.na(matched), NA, data[matched])
}

combined_df$DESeq2_pval <- safe_match(dat_list$DESeq2$pvalue, rownames(deseq2_res), all_genes)
combined_df$DESeq2_logFC <- safe_match(dat_list$DESeq2$logFC, rownames(deseq2_res), all_genes)
combined_df$eSVD_pval <- safe_match(dat_list$eSVD$pvalue, dat_list$eSVD$genes, all_genes)
combined_df$eSVD_logFC <- safe_match(dat_list$eSVD$logFC, dat_list$eSVD$genes, all_genes)
combined_df$NEBULA_pval <- safe_match(dat_list$NEBULA$pvalue, dat_list$NEBULA$genes, all_genes)
combined_df$NEBULA_logFC <- safe_match(dat_list$NEBULA$logFC, dat_list$NEBULA$genes, all_genes)
# combined_df$Was2_pval <- safe_match(dat_list$WAS2$pvalue, dat_list$WAS2$genes, all_genes)
# combined_df$Was2_logFC <- safe_match(dat_list$WAS2$logFC, dat_list$WAS2$genes, all_genes)
combined_df$test_var_pval <- safe_match(result_df$P_Value, result_df$Gene, all_genes)
combined_df$test_var_logFC <- safe_match(result_df$log2FC, result_df$Gene, all_genes)

output_path <- "~/kzlinlab/projects/subject-de/out/tati/Writeup7/ROSMAP_dataset_ingredients.csv"
write.csv(combined_df, file = output_path, row.names = FALSE)
