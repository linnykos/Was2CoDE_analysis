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

load("~/kzlinlab/projects/subject-de/out/tati/Writeup3/Writeup3_NEBULA.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup3/Writeup3_Katie_Pseudobulk-DEseq2.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup3/Writeup3_prater_esvd.RData")

dat_list <- list(
  eSVD = list(
    pvalue = 10^(-eSVD_obj$pvalue_list$log10pvalue),
    logFC = log2(eSVD_obj$case_mean / eSVD_obj$control),
    genes = names(eSVD_obj$teststat_vec)
  ),
  DESeq2 = list(
    pvalue = deseq2_res[,"pvalue"],
    logFC = deseq2_res[,"log2FoldChange"],
    genes = rownames(deseq2_res)
  ),
  NEBULA = list(
    pvalue = nebula_res$summary[,"p_Study_DesignationCtrl"],
    logFC = nebula_res$summary[,"logFC_Study_DesignationCtrl"],
    genes = nebula_res$summary$gene
  )
)

method_names <- names(dat_list)
combinations <- combn(method_names, 2, simplify = FALSE)

plot_combination <- function(comb) {
  method1 <- comb[1]
  method2 <- comb[2]
  data1 <- dat_list[[method1]]
  data2 <- dat_list[[method2]]
  
  intersect_genes <- intersect(data1$genes, data2$genes)
  indices1 <- match(intersect_genes, data1$genes)
  indices2 <- match(intersect_genes, data2$genes)
  
  # Fetch correct p-values and logFC using matched indices
  res <- data.frame(
    gene = intersect_genes,
    pvalue1 = data1$pvalue[indices1],
    logFC1 = data1$logFC[indices1],
    pvalue2 = data2$pvalue[indices2],
    logFC2 = data2$logFC[indices2]
  )
  ############################
  # Define thresholds
  ############################
  
  # Calculate adjusted p-values and define cutoffs
  res$pvalue_adj1 <- stats::p.adjust(res$pvalue1, method = "BH")
  res$pvalue_adj2 <- stats::p.adjust(res$pvalue2, method = "BH")
  idx1 <- which(res$pvalue_adj1 <= 0.05)
  idx2 <- which(res$pvalue_adj2 <= 0.05)
  pCutoff1 <- max(res[,"pvalue1"][idx1])
  pCutoff2 <- max(res[,"pvalue2"][idx2])
  
  # Compute thresholds for significant logFC
  FCcutoff1 <- quantile(abs(res$logFC1), 0.9, na.rm = TRUE)
  FCcutoff2 <- quantile(abs(res$logFC2), 0.9, na.rm = TRUE)
  
  # Set axis limits based on the 99th percentile of the absolute logFC values
  xlims <- c(-1, 1) * max(quantile(abs(res$logFC1), 0.99, na.rm = TRUE),
                          quantile(abs(res$logFC2), 0.99, na.rm = TRUE))

  ############################
  # Reorder the levels of SignificanceCategory
  ############################
  # Define significance categories based on adjusted p-values
  res$SignificanceCategory <- ifelse(-log10(res$pvalue1) > -log10(pCutoff1) & -log10(res$pvalue2) > -log10(pCutoff2), "Both",
                                     ifelse(-log10(res$pvalue1) > -log10(pCutoff1), method1,
                                            ifelse(-log10(res$pvalue2) > -log10(pCutoff2), method2, "Neither")))

  res$SignificanceCategory <- factor(res$SignificanceCategory, levels = c("Neither", method1, method2, "Both"))
  
  res <- res[c(which(res$SignificanceCategory == "Neither"),
               which(res$SignificanceCategory == method1),
               which(res$SignificanceCategory == method2),
               which(res$SignificanceCategory == "Both")),
  ]
  # vector for significant genes
  significant_genes <- res$SignificanceCategory != "Neither"
  correlation_pvalue <- stats::cor(-log10(res$pvalue1), -log10(res$pvalue2), use = "complete.obs")
  correlation_logfc <- stats::cor(res$logFC1, res$logFC2, use = "complete.obs")
  
  ############################
  # Generate Plots
  ############################
  total_genes <- nrow(res)
  significant_method1 <- sum(res$SignificanceCategory %in% c(method1, "Both"))
  significant_method2 <- sum(res$SignificanceCategory %in% c(method2, "Both"))
  significant_both <- sum(res$SignificanceCategory == "Both")
  
  colors <- setNames(c("#e0e0e0", "#cc0967", "#109163", "#bc6a17"), c("Neither", as.character(method1), as.character(method2), "Both"))
  
  plot_pvalue <- ggplot(res, aes(x = -log10(pvalue1), y = -log10(pvalue2), color = SignificanceCategory)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = colors) +
    geom_text_repel(aes(label = ifelse(SignificanceCategory != "Neither", as.character(gene), "")), box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = 10) +
    geom_vline(xintercept = -log10(pCutoff1), linetype = "dashed", color = "#cc0967") +
    geom_hline(yintercept = -log10(pCutoff2), linetype = "dashed", color = "#bc6a17") +
    labs(title = sprintf("P-value Comparison: %s vs. %s\nCorrelation: %.2f\nTotal Genes: %d, %s: %d, %s: %d, Both: %d",
                         method1, method2, stats::cor(-log10(res$pvalue1), -log10(res$pvalue2), use = "complete.obs"),
                         total_genes, method1, significant_method1, method2, significant_method2, significant_both),
         x = sprintf("-Log10 P-value (%s)", method1),
         y = sprintf("-Log10 P-value (%s)", method2)) +
    theme_minimal()
  
  plot_logfc <- ggplot(res, aes(x = logFC1, y = logFC2, color = SignificanceCategory)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = colors) +
    geom_text_repel(aes(label = ifelse(SignificanceCategory != "Neither", as.character(gene), "")), box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = 10) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
    geom_vline(xintercept = c(-FCcutoff1, FCcutoff1), linetype = "dotted", color = "#010086") +
    geom_hline(yintercept = c(-FCcutoff2, FCcutoff2), linetype = "dotted", color = "#bc6a17") +
    labs(title = sprintf("Log2FC Comparison: %s vs. %s\nCorrelation: %.2f", method1, method2, correlation_logfc),
         x = sprintf("Log2 FC (%s)", method1),
         y = sprintf("Log2 FC (%s)", method2)) +
    theme_minimal()
  
  # Save combined plots
  filename <- paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup3/Writeup3_Katie_Comparison_", method1, "to_", method2, ".png")
  combined_plot <- grid.arrange(plot_pvalue, plot_logfc, ncol = 2)
  ggsave(filename, combined_plot, device = "png", width = 14, height = 7, units = "in")
}

# Apply the function to each method combination
plots <- lapply(combinations, plot_combination)

