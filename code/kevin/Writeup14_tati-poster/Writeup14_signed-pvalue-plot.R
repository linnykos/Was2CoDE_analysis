rm(list=ls())

library(ggplot2)

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup14/"


.plot_ingredients <- function(df, method1, method2){
  x_lfc <- df[,paste0(method1, "_logFC")]
  x_sign <- sign(x_lfc)
  x_pval <- df[,paste0(method1, "_pval")]
  x_val <- x_sign * -log10(x_pval)
  names(x_val) <- df$Gene
  x_pval_adj <- stats::p.adjust(x_pval, method = "BH")
  x_fdrcutoff <- -log10(max(x_pval[which(x_pval_adj <= 0.05)]))
  if(is.nan(x_fdrcutoff)) x_fdrcutoff <- 0
  
  y_lfc <- df[,paste0(method2, "_logFC")]
  y_sign <- sign(y_lfc)
  y_pval <- df[,paste0(method2, "_pval")]
  y_val <- y_sign * -log10(y_pval)
  names(y_val) <- df$Gene
  y_pval_adj <- stats::p.adjust(y_pval, method = "BH")
  y_fdrcutoff <- -log10(max(y_pval[which(y_pval_adj <= 0.05)]))
  if(is.nan(y_fdrcutoff)) y_fdrcutoff <- 0
  
  # compute plotting ingredients
  idx <- intersect(which(!is.na(x_val)), which(!is.na(y_val)))
  x_val <- x_val[idx]; y_val <- y_val[idx]
  cor_val <- stats::cor(x_val, y_val)
  limit_vec <- c(-1,1)*max(abs(c(x_val, y_val)))
  
  # Perform PCA to get the leading principal component
  tmp_mat <- cbind(x_val, y_val)
  pca <- stats::prcomp(tmp_mat, center = FALSE)
  # Get the first principal component direction
  pc1_slope <- pca$rotation[2, 1] / pca$rotation[1, 1]
  
  # organize the plotting data frame
  ggplot_df <- data.frame(
    gene = names(x_val),
    method1 = x_val,
    method2 = y_val
  )
  
  # RGB color of purple: 
  purple_color <- rgb(135, 50, 255, maxColorValue = 255)
  
  # make the plot
  plot1 <- ggplot2::ggplot(data = ggplot_df, mapping = aes(x = method1,
                                                           y = method2)) + 
    geom_rect(data=ggplot_df[1,], aes(xmin = x_fdrcutoff, xmax = limit_vec[2], ymin = y_fdrcutoff, ymax = limit_vec[2]), fill = purple_color, alpha = 0.2) +
    geom_rect(data=ggplot_df[1,], aes(xmin = limit_vec[1], xmax = -x_fdrcutoff, ymin = limit_vec[1], ymax = -y_fdrcutoff), fill = purple_color, alpha = 0.2) +
    geom_point(color = "gray30", alpha = 0.5) +
    geom_hline(yintercept = y_fdrcutoff, linetype = "dashed") +  
    geom_hline(yintercept = -y_fdrcutoff, linetype = "dashed") +  
    geom_vline(xintercept = x_fdrcutoff, linetype = "dashed") + 
    geom_vline(xintercept = -x_fdrcutoff, linetype = "dashed") +
    geom_abline(slope = pc1_slope, intercept = 0, color = "coral", linewidth = 1) +
    xlim(limit_vec) + 
    ylim(limit_vec) +  
    scale_x_continuous(limits = limit_vec, expand = c(0, 0)) +  # Set exact x-limits without expansion
    scale_y_continuous(limits = limit_vec, expand = c(0, 0)) +  # Set exact y-limits without expansion
    coord_fixed(ratio = 1) +
    labs(x = paste0(method1, " (Signed -log10 p-value)"), 
         y = paste0(method2, " (Signed -log10 p-value)")) + 
    ggtitle(paste0(method1, " vs. ", method2, " (Cor: ", round(cor_val, 2), ")")) + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),        # Smaller title text size
      axis.title.x = element_text(size = 10),      # Smaller x-axis label text size
      axis.title.y = element_text(size = 10)       # Smaller y-axis label text size
    ) 
}

###############

df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup6/SEA-AD_dataset_ingredients.csv")

dataset_name <- "SEA-AD"
method1 <- "eSVD"
method2 <- "NEBULA"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

method1 <- "eSVD"
method2 <- "DESeq2"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

method1 <- "DESeq2"
method2 <- "NEBULA"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

#####################

df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Prater_dataset_ingredients.csv")

dataset_name <- "Prater"
method1 <- "eSVD"
method2 <- "NEBULA"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

method1 <- "eSVD"
method2 <- "DESeq2"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

method1 <- "DESeq2"
method2 <- "NEBULA"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

############


df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup7/ROSMAP_dataset_ingredients.csv")

dataset_name <- "ROSMAP"
method1 <- "eSVD"
method2 <- "NEBULA"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

method1 <- "eSVD"
method2 <- "DESeq2"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)

method1 <- "DESeq2"
method2 <- "NEBULA"
plot1 <- .plot_ingredients (df, 
                            method1 = method1, 
                            method2 = method2)
ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_", dataset_name, "_signed-pvalue_", method1, "-vs-", method2, ".png"),
       height = 3,
       width = 3)