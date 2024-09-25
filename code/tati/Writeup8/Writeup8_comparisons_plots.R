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

generate_comparison_plot <- function(data1, data2, method, dataset1_name, dataset2_name, output_dir) {
  
  # Find intersecting genes between the two datasets
  intersect_genes <- intersect(data1$genes, data2$genes)
  indices_data1 <- match(intersect_genes, data1$genes)
  indices_data2 <- match(intersect_genes, data2$genes)
  
  # Create a dataframe for plotting
  comparison_df <- data.frame(
    gene = intersect_genes,
    logFC_data1 = data1$logFC[indices_data1],
    logFC_data2 = data2$logFC[indices_data2],
    pvalue_data1 = data1$pvalue[indices_data1],
    pvalue_data2 = data2$pvalue[indices_data2]
  )
  
  # Check if any valid data is available for correlation and plotting
  if (nrow(comparison_df) == 0 || all(is.na(comparison_df$logFC_data1)) || all(is.na(comparison_df$logFC_data2))) {
    message("No valid logFC data found for method: ", method, " between ", dataset1_name, " and ", dataset2_name)
    return(NULL)  # Skip this plot if no valid data
  }
  
  # Filter out rows with missing logFC values
  comparison_df <- comparison_df[complete.cases(comparison_df$logFC_data1, comparison_df$logFC_data2), ]
  
  # If after filtering there is no data, return NULL
  if (nrow(comparison_df) == 0) {
    message("No complete cases after filtering for method: ", method, " between ", dataset1_name, " and ", dataset2_name)
    return(NULL)
  }
  
  # Define significance threshold for coloring with more granular categorization
  comparison_df$Significance <- ifelse(
    comparison_df$pvalue_data1 <= 0.05 & comparison_df$pvalue_data2 <= 0.05, 
    "Both Significant", 
    ifelse(
      comparison_df$pvalue_data1 <= 0.05, 
      sprintf("Only %s Significant", dataset1_name), 
      ifelse(
        comparison_df$pvalue_data2 <= 0.05, 
        sprintf("Only %s Significant", dataset2_name), 
        "Neither Significant"
      )
    )
  )
  
  # Calculate the correlation between logFC values, ignoring NA values
  correlation_logfc <- stats::cor(comparison_df$logFC_data1, comparison_df$logFC_data2, use = "complete.obs")
  
  # Define axis limits for the plot based on the 99th percentile of absolute logFC values
  xlims <- c(-1, 1) * max(quantile(abs(comparison_df$logFC_data1), 0.99, na.rm = TRUE),
                          quantile(abs(comparison_df$logFC_data2), 0.99, na.rm = TRUE))
  
  # Create the labels for significance categories
  sig_label_1 <- sprintf("Only %s Significant", dataset1_name)
  sig_label_2 <- sprintf("Only %s Significant", dataset2_name)
  
  # Generate the LogFC comparison plot
  plot_logfc <- ggplot(comparison_df, aes(x = logFC_data1, y = logFC_data2, color = Significance)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c(
      "Neither Significant" = "#2f2f2f", 
      sig_label_1 = "#109163",  # Explicit label for dataset1_name
      sig_label_2 = "#cc0967",  # Explicit label for dataset2_name
      "Both Significant" = "#bc6a17"
    )) +  
    geom_text_repel(aes(label = ifelse(Significance != "Neither Significant", as.character(gene), "")),
                    box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = 10) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
    labs(title = sprintf("Log2FC Comparison: %s on %s vs %s\nCorrelation: %.2f", 
                         method, dataset1_name, dataset2_name, correlation_logfc),
         x = sprintf("Log2 FC (%s)", dataset1_name), 
         y = sprintf("Log2 FC (%s)", dataset2_name)) +
    theme_minimal() +
    xlim(xlims[1], xlims[2]) +
    ylim(xlims[1], xlims[2])
  
  # Save plot
  filename <- sprintf("%s/Writeup8_Comparison_%s_%s_vs_%s_%s.png", output_dir, method, dataset1_name, dataset2_name, method)
  ggsave(filename, plot_logfc, device = "png", width = 7, height = 7, units = "in")
  
  # Return the plot object (optional)
  return(plot_logfc)
}

######################################################

# Dataset lists for each method
load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_NEBULA.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_Pseudobulk-DEseq2.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_esvd.RData")

dat_list_prater <- list(
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
    pvalue = nebula_res$summary[,"p_Study_DesignationAD"],
    logFC = nebula_res$summary[,"logFC_Study_DesignationAD"],
    genes = nebula_res$summary$gene
  )
)
load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_NEBULA.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_Pseudobulk-DEseq2.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_eSVD.RData")

dat_list_seaad <- list(
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
    pvalue = nebula_res$summary[,"p_ADNCCase"],
    logFC = nebula_res$summary[,"logFC_ADNCCase"],
    genes = nebula_res$summary$gene
  )
)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_NEBULA.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_Pseudobulk-DEseq2.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_eSVD.RData")

dat_list_rosmap <- list(
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

# Define dataset pairs to compare
dataset_pairs <- list(
  list("Prater", "SEA-AD", dat_list_prater, dat_list_seaad),
  list("SEA-AD", "ROSMAP", dat_list_seaad, dat_list_rosmap),
  list("Prater", "ROSMAP", dat_list_prater, dat_list_rosmap)
)

# Define methods to compare
methods <- c("eSVD", "DESeq2", "NEBULA")

# Output directory for saving plots
output_dir <- "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8"

# Generate comparison plots for all methods and dataset pairs
for (pair in dataset_pairs) {
  dataset1_name <- pair[[1]]
  dataset2_name <- pair[[2]]
  data1 <- pair[[3]]
  data2 <- pair[[4]]
  
  for (method in methods) {
    plot <- generate_comparison_plot(data1[[method]], data2[[method]], method, dataset1_name, dataset2_name, output_dir)
  }
}

##########################################
# load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_NEBULA.RData")
# nebula_res_prater <- nebula_res
# load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_NEBULA.RData")
# nebula_res_SEAAD <- nebula_res
# # Extract logFC and gene names from prater's data (Writeup5) and SEA-AD data (Writeup6)
# prater_data <- list(
#   logFC = nebula_res_prater$summary$logFC_Study_DesignationAD,
#   genes = nebula_res_prater$summary$gene
# )
# 
# seaad_data <- list(
#   logFC = nebula_res_SEAAD$summary$logFC_ADNCCase,
#   genes = nebula_res_SEAAD$summary$gene
# )
# 
# # Find intersecting genes between the two datasets
# intersect_genes <- intersect(prater_data$genes, seaad_data$genes)
# indices_prater <- match(intersect_genes, prater_data$genes)
# indices_seaad <- match(intersect_genes, seaad_data$genes)
# 
# # Create a dataframe for plotting
# comparison_df <- data.frame(
#   gene = intersect_genes,
#   logFC_prater = prater_data$logFC[indices_prater],
#   logFC_seaad = seaad_data$logFC[indices_seaad]
# )
# 
# # Define significance thresholds for coloring (|logFC| > 1)
# comparison_df$Significance <- ifelse(abs(comparison_df$logFC_prater) > 1 | abs(comparison_df$logFC_seaad) > 1, "Significant", "Not Significant")
# 
# # Calculate the correlation between logFC values
# correlation_logfc <- stats::cor(comparison_df$logFC_prater, comparison_df$logFC_seaad, use = "complete.obs")
# 
# # Define axis limits for the plot based on the 99th percentile of absolute logFC values
# xlims <- c(-1, 1) * max(quantile(abs(comparison_df$logFC_prater), 0.99, na.rm = TRUE),
#                         quantile(abs(comparison_df$logFC_seaad), 0.99, na.rm = TRUE))
# 
# # Generate the LogFC comparison plot
# plot_logfc <- ggplot(comparison_df, aes(x = logFC_prater, y = logFC_seaad, color = Significance)) +
#   geom_point(alpha = 0.7) +
#   scale_color_manual(values = c("Not Significant" = "#2f2f2f", "Significant" = "#cc0967")) +  # Color points based on significance
#   geom_text_repel(aes(label = ifelse(Significance == "Significant", as.character(gene), "")),
#                   box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = 10) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
#   labs(title = sprintf("Log2FC Comparison: NEBULA on prater's vs SEA-AD\nCorrelation: %.2f", correlation_logfc),
#        x = "Log2 FC (prater's Data)", y = "Log2 FC (SEA-AD Data)") +
#   theme_minimal() +
#   xlim(xlims[1], xlims[2]) +
#   ylim(xlims[1], xlims[2])
# # Save plot
# filename <- "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8/Writeup8_Comparison_NEBULA_prater_vs_SEAAD.png"
# ggsave(filename, plot_logfc, device = "png", width = 7, height = 7, units = "in")
# 
# ############################
# 
# 
# 
# load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_Pseudobulk-DEseq2.RData")
# load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_eSVD.RData")
# load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_was2_wilcox.RData")
# 
# 
# dat_list <- list(
#   eSVD = list(
#     pvalue = 10^(-eSVD_obj$pvalue_list$log10pvalue),
#     logFC = log2(eSVD_obj$case_mean / eSVD_obj$control),
#     genes = names(eSVD_obj$teststat_vec)
#   ),
#   WAS2 = list(
#     pvalue = results_mat[,"p_val"],
#     logFC = log2((results_mat[, "mean_dn"]) / ((results_mat[, "mean_nn"] + results_mat[, "mean_dd"]) / 2)),
#     genes = rownames(results_mat)
#   ),
#   DESeq2 = list(
#     pvalue = deseq2_res[,"pvalue"],
#     logFC = deseq2_res[,"log2FoldChange"],
#     genes = rownames(deseq2_res)
#   ),
#   NEBULA = list(
#     pvalue = nebula_res$summary[,"p_ADNCCase"],
#     logFC = nebula_res$summary[,"logFC_ADNCCase"],
#     genes = nebula_res$summary$gene
#   )
# )
# 
# method_names <- names(dat_list)
# combinations <- combn(method_names, 2, simplify = FALSE)
# 
# plot_combination <- function(comb) {
#   method1 <- comb[1]
#   method2 <- comb[2]
#   data1 <- dat_list[[method1]]
#   data2 <- dat_list[[method2]]
#   
#   intersect_genes <- intersect(data1$genes, data2$genes)
#   indices1 <- match(intersect_genes, data1$genes)
#   indices2 <- match(intersect_genes, data2$genes)
#   
#   # Fetch correct p-values and logFC using matched indices
#   res <- data.frame(
#     gene = intersect_genes,
#     pvalue1 = data1$pvalue[indices1],
#     logFC1 = data1$logFC[indices1],
#     pvalue2 = data2$pvalue[indices2],
#     logFC2 = data2$logFC[indices2]
#   )
#   ############################
#   # Define thresholds
#   ############################
#   
#   # Calculate adjusted p-values and define cutoffs
#   res$pvalue_adj1 <- stats::p.adjust(res$pvalue1, method = "BH")
#   res$pvalue_adj2 <- stats::p.adjust(res$pvalue2, method = "BH")
#   idx1 <- which(res$pvalue_adj1 <= 0.05)
#   idx2 <- which(res$pvalue_adj2 <= 0.05)
#   pCutoff1 <- max(res[,"pvalue1"][idx1])
#   pCutoff2 <- max(res[,"pvalue2"][idx2])
#   
#   # Compute thresholds for significant logFC
#   FCcutoff1 <- quantile(abs(res$logFC1), 0.9, na.rm = TRUE)
#   FCcutoff2 <- quantile(abs(res$logFC2), 0.9, na.rm = TRUE)
#   
#   # Set axis limits based on the 99th percentile of the absolute logFC values
#   xlims <- c(-1, 1) * max(quantile(abs(res$logFC1), 0.99, na.rm = TRUE),
#                           quantile(abs(res$logFC2), 0.99, na.rm = TRUE))
#   
#   ############################
#   # Reorder the levels of SignificanceCategory
#   ############################
#   # Define significance categories based on adjusted p-values
#   res$SignificanceCategory <- ifelse(-log10(res$pvalue1) > -log10(pCutoff1) & -log10(res$pvalue2) > -log10(pCutoff2), "Both",
#                                      ifelse(-log10(res$pvalue1) > -log10(pCutoff1), method1,
#                                             ifelse(-log10(res$pvalue2) > -log10(pCutoff2), method2, "Neither")))
#   
#   res$SignificanceCategory <- factor(res$SignificanceCategory, levels = c("Neither", method1, method2, "Both"))
#   
#   res <- res[c(which(res$SignificanceCategory == "Neither"),
#                which(res$SignificanceCategory == method1),
#                which(res$SignificanceCategory == method2),
#                which(res$SignificanceCategory == "Both")),
#   ]
#   # vector for significant genes
#   significant_genes <- res$SignificanceCategory != "Neither"
#   correlation_pvalue <- stats::cor(-log10(res$pvalue1), -log10(res$pvalue2), use = "complete.obs")
#   correlation_logfc <- stats::cor(res$logFC1, res$logFC2, use = "complete.obs")
#   
#   ############################
#   # Generate Plots
#   ############################
#   total_genes <- nrow(res)
#   significant_method1 <- sum(res$SignificanceCategory %in% c(method1, "Both"))
#   significant_method2 <- sum(res$SignificanceCategory %in% c(method2, "Both"))
#   significant_both <- sum(res$SignificanceCategory == "Both")
#   
#   colors <- setNames(c("#e0e0e0", "#cc0967", "#109163", "#bc6a17"), c("Neither", as.character(method1), as.character(method2), "Both"))
#   
#   plot_pvalue <- ggplot(res, aes(x = -log10(pvalue1), y = -log10(pvalue2), color = SignificanceCategory)) +
#     geom_point(alpha = 0.7) +
#     scale_color_manual(values = colors) +
#     geom_text_repel(aes(label = ifelse(SignificanceCategory != "Neither", as.character(gene), "")), box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = 10) +
#     geom_vline(xintercept = -log10(pCutoff1), linetype = "dashed", color = "#cc0967") +
#     geom_hline(yintercept = -log10(pCutoff2), linetype = "dashed", color = "#bc6a17") +
#     labs(title = sprintf("P-value Comparison: %s vs. %s\nCorrelation: %.2f\nTotal Genes: %d, %s: %d, %s: %d, Both: %d",
#                          method1, method2, stats::cor(-log10(res$pvalue1), -log10(res$pvalue2), use = "complete.obs"),
#                          total_genes, method1, significant_method1, method2, significant_method2, significant_both),
#          x = sprintf("-Log10 P-value (%s)", method1),
#          y = sprintf("-Log10 P-value (%s)", method2)) +
#     theme_minimal()
#   
#   plot_logfc <- ggplot(res, aes(x = logFC1, y = logFC2, color = SignificanceCategory)) +
#     geom_point(alpha = 0.7) +
#     scale_color_manual(values = colors) +
#     geom_text_repel(aes(label = ifelse(SignificanceCategory != "Neither", as.character(gene), "")), box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = 10) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#2f2f2f") +
#     geom_vline(xintercept = c(-FCcutoff1, FCcutoff1), linetype = "dotted", color = "#010086") +
#     geom_hline(yintercept = c(-FCcutoff2, FCcutoff2), linetype = "dotted", color = "#bc6a17") +
#     labs(title = sprintf("Log2FC Comparison: %s vs. %s\nCorrelation: %.2f", method1, method2, correlation_logfc),
#          x = sprintf("Log2 FC (%s)", method1),
#          y = sprintf("Log2 FC (%s)", method2)) +
#     theme_minimal()
#   
#   # Save combined plots
#   filename <- paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup6/Writeup6_SEA-AD_Comparison_", method1, "_to_", method2, ".png")
#   combined_plot <- grid.arrange(plot_pvalue, plot_logfc, ncol = 2)
#   ggsave(filename, combined_plot, device = "png", width = 14, height = 7, units = "in")
# }
# 
# # Apply the function to each method combination
# plots <- lapply(combinations, plot_combination)
# 
