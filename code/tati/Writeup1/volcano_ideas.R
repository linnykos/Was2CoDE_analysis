rm(list=ls())

library(EnhancedVolcano)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideas.RData")

#compute avg_log2FC
dementia_indices <- which(meta_ind$CognitiveStatus == "Dementia")
no_dementia_indices <- which(meta_ind$CognitiveStatus == "No_dementia")
# create a vector to store avg_log2FC 
avg_log2FC <- numeric(dim(dist1)[1])
for (i in 1:dim(dist1)[1]) {
  gene_submatrix <- dist1[i, , ]
  avg_expr_dementia <- mean(gene_submatrix[dementia_indices, dementia_indices], na.rm = TRUE)
  avg_expr_no_dimentia <- mean(gene_submatrix[no_dementia_indices, no_dementia_indices], na.rm = TRUE)
  # Compute the average log fold change
  avg_log2FC[i] <- log2(avg_expr_dementia / avg_expr_no_dimentia)
}

quantile(pval_ideas)

volcano_data <- data.frame(
  gene = rownames(dist1),  
  avg_log2FC,
  pval_ideas)
# determine p-value cutoff
tmp <- sort(pval_ideas, decreasing = F)
tmp <- tmp[tmp != 0]
pCutoff <- tmp[10]

# determine log fold change cutoff
FCcutoff <- sort(abs(avg_log2FC), decreasing = TRUE)[100]

plot1 <- EnhancedVolcano::EnhancedVolcano(volcano_data,
                         lab = rownames(volcano_data),
                         x = "avg_log2FC",
                         y = "pval_ideas",
                         pCutoff = pCutoff,
                         FCcutoff = FCcutoff,
                         title = 'Volcano Plot',
                         subtitle = 'Comparison of Dementia vs. No Dementia',
                         xlab = 'Log2 Fold Change',
                         ylab = '-Log10 pval_ideas',
                         #selectLab = volcano_data$gene[volcano_data$pval_ideas < pCutoff & abs(volcano_data$avg_log2FC) > FCcutoff],
                         #selectLab = volcano_data$gene[which((volcano_data$pval_ideas < pCutoff) | (abs(volcano_data$avg_log2FC) > FCcutoff)])
                         )

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-volcano_microglia_ideas_Lab.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

###########


# determine p-value cutoff
tmp <- sort(de_results_t$p_val_adj, decreasing = F)
idx <- which.max(de_results_t$p_val_adj[de_results_t$p_val_adj <= 0.05])
pCutoff <- de_results_t$p_val[idx]

# determine log fold change cutoff
FCcutoff <- sort(abs(de_results_t$avg_log2FC), decreasing = T)[100]

plot1 <-  EnhancedVolcano::EnhancedVolcano(de_results_t,
                                           lab = rownames(de_results_t),
                                           x = "avg_log2FC",
                                           y = "p_val_adj",
                                           pCutoff = pCutoff,
                                           FCcutoff = FCcutoff)
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-volcano_t.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")