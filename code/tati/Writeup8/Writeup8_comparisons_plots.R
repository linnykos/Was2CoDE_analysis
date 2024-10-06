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
generate_comparison_plot2 <- function(data1, data2, method, dataset1_name, dataset2_name, output_dir) {
  intersect_genes <- intersect(data1$genes, data2$genes)
  indices_data1 <- match(intersect_genes, data1$genes)
  indices_data2 <- match(intersect_genes, data2$genes)
  
  # Create a dataframe for plotting
  comparison_df <- data.frame(
    gene = intersect_genes,
    log10p_data1 = -log10(data1$pvalue[indices_data1]),
    log10p_data2 = -log10(data2$pvalue[indices_data2]),
    logFC_data1 = data1$logFC[indices_data1],
    logFC_data2 = data2$logFC[indices_data2]
  )
  
  # Remove rows with infinite values
  comparison_df <- comparison_df[is.finite(comparison_df$log10p_data1) & 
                                   is.finite(comparison_df$log10p_data2), ]
  
  # Define significance thresholds (you may need to adjust these)
  thresh_x <- -log10(0.05)
  thresh_y <- -log10(0.05)
  
  # Determine point colors based on quadrants
  comparison_df$color <- case_when(
    comparison_df$log10p_data1 > thresh_x & comparison_df$log10p_data2 > thresh_y ~ "red",
    comparison_df$log10p_data1 > thresh_x & comparison_df$log10p_data2 <= thresh_y ~ "orange",
    comparison_df$log10p_data1 <= thresh_x & comparison_df$log10p_data2 > thresh_y ~ "blue",
    TRUE ~ "gray"
  )
  
  # correlation
  correlation <- cor(comparison_df$log10p_data1, comparison_df$log10p_data2, use = "complete.obs")
  
  p <- ggplot(comparison_df, aes(x = log10p_data1, y = log10p_data2)) +
    geom_point(aes(color = color), alpha = 0.7) +
    scale_color_identity() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_vline(xintercept = thresh_x, linetype = "dashed") +
    geom_hline(yintercept = thresh_y, linetype = "dashed") +
    labs(
      title = paste0(method, ": ", dataset1_name, " vs ", dataset2_name),
      x = paste0("-log10(p) (", dataset1_name, ")"),
      y = paste0("-log10(p) (", dataset2_name, ")")
    ) +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(0, max(comparison_df$log10p_data1, comparison_df$log10p_data2)),
                ylim = c(0, max(comparison_df$log10p_data1, comparison_df$log10p_data2)))
  
  # Add gene labels for top genes
  top_genes <- comparison_df %>%
    filter(log10p_data1 > thresh_x | log10p_data2 > thresh_y) %>%
    arrange(desc(pmax(log10p_data1, log10p_data2))) %>%
    head(10)
  
  p <- p + geom_text_repel(
    data = top_genes,
    aes(label = gene),
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = 'grey50'
  )
  
  # Add correlation to the plot
  p <- p + annotate("text", x = Inf, y = -Inf, 
                    label = sprintf("Correlation: %.2f", correlation),
                    hjust = 1, vjust = 0, size = 3)
  
  # Save the plot
  filename <- file.path(output_dir, sprintf("%s_%s_vs_%s_log10p_comparison.png", method, dataset1_name, dataset2_name))
  ggsave(filename, p, width = 10, height = 10, dpi = 300)
  
  return(p)
}

# Generate comparison plots for all methods and dataset pairs
for (pair in dataset_pairs) {
  dataset1_name <- pair[[1]]
  dataset2_name <- pair[[2]]
  data1 <- pair[[3]]
  data2 <- pair[[4]]
  
  for (method in methods) {
    plot <- generate_comparison_plot2(data1[[method]], data2[[method]], method, dataset1_name, dataset2_name, output_dir)
  }
}
##########################################
