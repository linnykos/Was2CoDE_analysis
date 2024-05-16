rm(list=ls())

library(EnhancedVolcano)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideascustom.RData")


#compute avg_log2FC
dementia_indices <- which(meta_ind$CognitiveStatus == "Dementia")
no_dementia_indices <- which(meta_ind$CognitiveStatus == "No_dementia")
# create a vector to store avg_log2FC 
for (i in 1:6){
  avg_log2FC <- numeric(dim(dist_list[[i]])[1])
  
  for (k in 1:dim(dist_list[[i]])[1]) {
    gene_submatrix <- dist_list[[i]][k, , ]
    avg_expr_dementia <- mean(gene_submatrix[dementia_indices, dementia_indices], na.rm = TRUE)
    avg_expr_no_dimentia <- mean(gene_submatrix[no_dementia_indices, no_dementia_indices], na.rm = TRUE)
    
    # Compute the average log fold change
    avg_log2FC[k] <- log2(avg_expr_dementia / avg_expr_no_dimentia)
  }

  quantile(pval_ideas[[i]])
  volcano_data <- data.frame(
    gene = rownames(count_matrix),  
    avg_log2FC,
    pval = pval_ideas[[i]])
  
  
  
  # determine p-value cutoff
  tmp <- sort(pval_ideas[[i]], decreasing = F)
  tmp <- tmp[tmp != 0]
  pCutoff <- tmp[10]
  
  # determine log fold change cutoff
  FCcutoff <- sort(abs(avg_log2FC), decreasing = TRUE)[100]
  Plot1<- EnhancedVolcano::EnhancedVolcano(volcano_data,
                                                lab = rownames(volcano_data),
                                                x = "avg_log2FC",
                                                y = "pval",
                                                pCutoff = pCutoff,
                                                FCcutoff = FCcutoff,
                                                title = 'Volcano Plot',
                                                subtitle = 'Comparison of Dementia vs. No Dementia',
                                                xlab = 'Log2 Fold Change',
                                                ylab = '-Log10 pval_ideas')
                                          
}


ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-volcano_microglia_IdeasCustom_1_Labeled2.png"),
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
###

dementia_indices <- which(meta_ind$CognitiveStatus == "Dementia")
no_dementia_indices <- which(meta_ind$CognitiveStatus == "No_dementia")
# create a vector to store avg_log2FC 
avg_log2FC <- numeric(dim(dist_list[[1]])[1])

for (k in 1:dim(dist_list[[1]])[1]) {
gene_submatrix <- dist_list[[1]][k, , ]
avg_expr_dementia <- mean(gene_submatrix[dementia_indices, dementia_indices], na.rm = TRUE)
avg_expr_no_dimentia <- mean(gene_submatrix[no_dementia_indices, no_dementia_indices], na.rm = TRUE)

 # Compute the average log fold change
 avg_log2FC[k] <- log2(avg_expr_dementia / avg_expr_no_dimentia)
}
quantile(pval_ideas[[1]])
volcano_data <- data.frame(
  gene = rownames(count_matrix),  
  avg_log2FC,
  pval = pval_ideas[[1]])
# determine p-value cutoff
tmp <- sort(pval_ideas[[1]], decreasing = F)
tmp <- tmp[tmp != 0]
pCutoff <- tmp[10]
FCcutoff <- sort(abs(avg_log2FC), decreasing = TRUE)[100]
volcano_data$pval <- as.numeric(volcano_data$pval)
plot1 <- EnhancedVolcano::EnhancedVolcano(volcano_data,
                                          lab = volcano_data$gene,
                                          x = "avg_log2FC",
                                          y = "pval",
                                          pCutoff = pCutoff,
                                          FCcutoff = FCcutoff,
                                          title = 'Volcano Plot',
                                          subtitle = 'Comparison of Dementia vs. No Dementia',
                                          xlab = 'Log2 Fold Change',
                                          ylab = '-Log10 pval_ideas')
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-volcano_microglia_IdeasCustom_1_Labeled4.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")
plot1 <- EnhancedVolcano::EnhancedVolcano(volcano_data,
                                          lab = rownames(volcano_data),
                                          x = "avg_log2FC",
                                          y = "pval",
                                          pCutoff = pCutoff,
                                          FCcutoff = FCcutoff,
                                          title = 'Volcano Plot',
                                          subtitle = 'Comparison of Dementia vs. No Dementia',
                                          xlab = 'Log2 Fold Change',
                                          ylab = '-Log10 pval_ideas')
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-volcano_microglia_IdeasCustom_1_Labeled5.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")