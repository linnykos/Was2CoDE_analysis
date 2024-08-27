rm(list=ls())
library(ggplot2)
library(ggrepel)
library(dplyr)
library(Seurat)

set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_ideascustom.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_Katie_was2_wilcox.RData")
ls()

compute_statistics <- function(dist_list, meta_ind, pval_list) {
  n_gene <- dim(dist_list$distance)[1]
  n_donors <- dim(dist_list$distance)[2]
  
  statistics <- matrix(NA, nrow = n_gene, ncol = 5)
  colnames(statistics) <- c("median_loc_dist2_case_control", 
                            "median_loc_dist2_case_case_control_control", 
                            "median_size_dist2_case_control", 
                            "median_size_dist2_case_case_control_control", 
                            "-log10_p_value")
  rownames(statistics) <- names(pval_list)
  
  for (i in 1:n_gene) {
    location <- dist_list$location[i, , ]
    size <- dist_list$size[i, , ]
    distance <- dist_list$distance[i, , ]
    
    loc_dist2 <- location / (distance^2)
    size_dist2 <- size / (distance^2)
    
    case_indices <- which(meta_ind$CognitiveStatus == "Dementia")
    control_indices <- which(meta_ind$CognitiveStatus == "No_dementia")
    
    loc_dist2_case_control <- loc_dist2[case_indices, control_indices]
    size_dist2_case_control <- size_dist2[case_indices, control_indices]
    
    loc_dist2_case_case_control_control <- loc_dist2[c(case_indices, control_indices), c(case_indices, control_indices)]
    size_dist2_case_case_control_control <- size_dist2[c(case_indices, control_indices), c(case_indices, control_indices)]
    
    statistics[i, 1] <- median(loc_dist2_case_control, na.rm = TRUE)
    statistics[i, 2] <- median(loc_dist2_case_case_control_control, na.rm = TRUE)
    statistics[i, 3] <- median(size_dist2_case_control, na.rm = TRUE)
    statistics[i, 4] <- median(size_dist2_case_case_control_control, na.rm = TRUE)
    
    # Assuming we want to use the p-value for the 'distance' metric for the log transformation
    p_value <- pval_list[i]
    # names(p_value) <- names(pval_list[i])
    statistics[i, 5] <- -log10(p_value)
    # names(statistics[i, 5]) <- names(p_value)
  }
  
  return(statistics)
}

statistics <- compute_statistics(dist_list, meta_ind, p_values)

statistics_df <- as.data.frame(statistics)
# Arrange the data frame by -log10_p_value in descending order
colnames(statistics_df) <- c("median_loc_dist2_case_control", 
                             "median_loc_dist2_case_case_control_control", 
                             "median_size_dist2_case_control", 
                             "median_size_dist2_case_case_control_control", 
                             "log10pvalue")
# Arrange the data frame by -log10_p_value in descending order
statistics_df <- statistics_df %>%
  arrange(`log10pvalue`)
# Labels
statistics_df$label <- ifelse(rank(-statistics_df$"log10pvalue") <= 10, rownames(statistics_df), "")
# To get the 10 genes with highest pvals
# tail(statistics_df[,5:6],10)
mid <- mean(statistics_df$`log10pvalue`, na.rm = TRUE)

# Plot 1: Median of location/distance^2 for pairs of donors case-control vs
# Median of location/distance^2 for pairs of donors case-case or control-control
# Plot 1: Med location/distance^2 (Case-Control vs Case-Case or Control-Control)
plot1 <- ggplot(statistics_df, aes(x = median_loc_dist2_case_case_control_control, 
                                   y = median_loc_dist2_case_control, 
                                   label = label, 
                                   color = `log10pvalue`)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(max.overlaps = 200,, size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 1: Med location/distance^2 (Case-Control vs Case-Case or Control-Control)",
       x = "Med location/dist^2 (Case-Case or Control-Control)",
       y = "Med location/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))
plot1 <-plot1+scale_color_gradient2(midpoint = mid, low = "blue", mid = "beige", high = "red", space = "Lab") 
# print(plot1)
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_was2_location_distance.png"),
                plot = plot1, device = "png", width = 8, height = 7, units = "in")


# Plot 2: Median of size/distance^2 (Case-Control vs Case-Case or Control-Control)
plot2 <- ggplot(statistics_df, aes(x = median_size_dist2_case_case_control_control, 
                                   y = median_size_dist2_case_control, 
                                   label = label, 
                                   color = `log10pvalue`)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(max.overlaps = 200) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 2: Median of size/distance^2 (Case-Control vs Case-Case or Control-Control)",
       x = "Median size/dist^2 (Case-Case or Control-Control)",
       y = "Median size/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))

plot2 <-plot2+scale_color_gradient2(midpoint = mid, low = "blue", mid = "beige", high = "red", space = "Lab")

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_was2_size_distance.png"),
                plot = plot2, device = "png", width = 8, height = 7, units = "in")

# Plot 3: Median of size/distance^2 vs location/dist^2 Case-Control
plot3 <- ggplot(statistics_df,
                aes(x = median_loc_dist2_case_control, 
                    y = median_size_dist2_case_control, 
                    label = label,
                    color = `log10pvalue`)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(max.overlaps = 200) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 3: size/distance^2 vs location/dist^2 Case-Control",
       x = "Median location/dist^2 (Case-Control)",
       y = "Median size/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "beige", high = "red", space = "Lab")

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_was2_size_location.png"),
                plot = plot3, device = "png", width = 8, height = 7, units = "in")

############ Volcano plot ############

library(EnhancedVolcano)

#compute avg_log2FC
dementia_indices <- which(meta_ind$CognitiveStatus == "Dementia")
no_dementia_indices <- which(meta_ind$CognitiveStatus == "No_dementia")
# create a vector to store avg_log2FC 
i=1
  avg_log2FC <- numeric(dim(dist_list[[i]])[1])
  
  for (k in 1:dim(dist_list[[i]])[1]) {
    gene_submatrix <- dist_list[[i]][k, , ]
    avg_expr_dementia <- mean(gene_submatrix[dementia_indices, dementia_indices], na.rm = TRUE)
    avg_expr_no_dimentia <- mean(gene_submatrix[no_dementia_indices, no_dementia_indices], na.rm = TRUE)
    
    # Compute the average log fold change
    avg_log2FC[k] <- log2(avg_expr_dementia / avg_expr_no_dimentia)
  }
  
  quantile(p_values)
  volcano_data <- data.frame(
    gene = names(p_values),  
    avg_log2FC,
    pval = p_values)
  
  # determine p-value cutoff
  tmp <- sort(p_values, decreasing = F)
  tmp <- tmp[tmp != 0]
  pCutoff <- tmp[10]
  
  # determine log fold change cutoff
  FCcutoff <- sort(abs(avg_log2FC), decreasing = TRUE)[100]
  Plot4<- EnhancedVolcano::EnhancedVolcano(volcano_data,
                                           lab = rownames(volcano_data),
                                           x = "avg_log2FC",
                                           y = "pval",
                                           pCutoff = pCutoff,
                                           FCcutoff = FCcutoff,
                                           title = 'Volcano Plot',
                                           subtitle = 'Comparison of Dementia vs. No Dementia',
                                           xlab = 'Log2 Fold Change',
                                           ylab = '-Log10 pval_ideas')



ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_was2_volcano.png"),
                Plot4, device = "png", width = 7, height = 7, units = "in")