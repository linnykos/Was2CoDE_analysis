rm(list=ls())
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

# from https://github.com/linnykos/subject-de/blob/tati/code/tati/Writeup1/Writeup1_microglia_ideascustom_5statistics.R
load("~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideascustom_5statistics.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup1/Writeup1_nebula.RData")

# Convert the statistics matrix to a data frame for ggplot2
statistics_df <- as.data.frame(statistics)
# Arrange the data frame by -log10_p_value in descending order
colnames(statistics_df) <- c("median_loc_dist2_case_control", 
                             "median_loc_dist2_case_case_control_control", 
                             "median_size_dist2_case_control", 
                             "median_size_dist2_case_case_control_control", 
                             "log10pvalue")

# Get NEBULA p-values
res <- nebula_res$summary
pval_vec <- res[,"p_CognitiveStatusNo dementia"]
names(pval_vec) <- res$gene
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
nebula_genes <- names(pval_adj_vec)[which(pval_adj_vec <= 0.05)]
nebula_genes <- intersect(nebula_genes, rownames(statistics_df))

# Labels
statistics_df$label <- rep("", nrow(statistics_df))
for(gene in nebula_genes){
  idx <- which(rownames(statistics_df) == gene)
  statistics_df$label[idx] <- gene
}
statistics_df$coloring <- rep(FALSE, nrow(statistics_df))
statistics_df$coloring[which(statistics_df$label != "")] <- TRUE

statistics_df <- statistics_df[c(which(statistics_df$label == ""),
                                 which(statistics_df$label != "")),]
# To get the 10 genes with highest pvals
mid <- mean(statistics_df$`log10pvalue`, na.rm = TRUE)

# Plot 1: Median of location/distance^2 for pairs of donors case-control vs
# Median of location/distance^2 for pairs of donors case-case or control-control
# Plot 1: Med location/distance^2 (Case-Control vs Case-Case or Control-Control)
plot1 <- ggplot(statistics_df, aes(x = median_loc_dist2_case_case_control_control, 
                                   y = median_loc_dist2_case_control, 
                                   label = label, 
                                   color = coloring)) 
plot1 <- plot1 + geom_point(alpha = 0.8) 
plot1 <- plot1 + geom_text_repel(max.overlaps = 200, 
                                 box.padding = ggplot2::unit(0.5, 'lines'),
                                 point.padding = ggplot2::unit(1.6, 'lines')) 
plot1 <- plot1 + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") 
plot1 <- plot1 + labs(title = "Plot 1: Med location/distance^2 (Case-Control vs Case-Case or Control-Control)",
       x = "Med location/dist^2 (Case-Case or Control-Control)",
       y = "Med location/dist^2 (Case-Control)",
       color = "-log10(p-value)") 
plot1 <- plot1 + theme_minimal() 
plot1 <- plot1 + theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
plot1 <- plot1 + xlim(c(0, 0.5)) + ylim(c(0, 0.5))
# print(plot1)
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_de_location_distance_wNebula.png"),
                plot = plot1, device = "png", width = 8, height = 7, units = "in")

# Plot 2: Median of size/distance^2 (Case-Control vs Case-Case or Control-Control)
plot2 <- ggplot(statistics_df, aes(x = median_size_dist2_case_case_control_control, 
                                   y = median_size_dist2_case_control, 
                                   label = label, 
                                   color = coloring)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(max.overlaps = 200,
                  box.padding = ggplot2::unit(0.5, 'lines'),
                  point.padding = ggplot2::unit(1.6, 'lines')) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 2: Median of size/distance^2 (Case-Control vs Case-Case or Control-Control)",
       x = "Median size/dist^2 (Case-Case or Control-Control)",
       y = "Median size/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))
plot2 <- plot2 + ggplot2::scale_colour_manual(values=c("black", "red"))
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_de_size_distance_wNebula.png"),
                plot = plot2, device = "png", width = 8, height = 7, units = "in")

# Plot 3: Median of size/distance^2 vs location/dist^2 Case-Control
plot3 <- ggplot(statistics_df,
                aes(x = median_loc_dist2_case_control, 
                    y = median_size_dist2_case_control, 
                    label = label,
                    color = coloring)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(max.overlaps = 200,
                  box.padding = ggplot2::unit(0.5, 'lines'),
                  point.padding = ggplot2::unit(1.6, 'lines')) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 3: size/distance^2 vs location/dist^2 Case-Control",
       x = "Median location/dist^2 (Case-Control)",
       y = "Median size/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) 
plot3 <- plot3 + ggplot2::scale_colour_manual(values=c("black", "red"))
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_de_size_location_wNebula.png"),
                plot = plot3, device = "png", width = 8, height = 7, units = "in")