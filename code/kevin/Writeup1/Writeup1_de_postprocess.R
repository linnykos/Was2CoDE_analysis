rm(list=ls())

library(EnhancedVolcano)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup1/Writeup1_de.RData")

dim(de_results_wilcoxon)
head(de_results_wilcoxon)
quantile(de_results_wilcoxon$p_val)

# determine p-value cutoff
tmp <- sort(de_results_wilcoxon$p_val, decreasing = F)
tmp <- tmp[tmp != 0]
pCutoff <- tmp[10]

# determine log fold change cutoff
FCcutoff <- sort(abs(de_results_wilcoxon$avg_log2FC), decreasing = T)[100]

plot1 <-  EnhancedVolcano::EnhancedVolcano(de_results_wilcoxon,
                                        lab = rownames(de_results_wilcoxon),
                                        x = "avg_log2FC",
                                        y = "p_val_adj",
                                        pCutoff = pCutoff,
                                        FCcutoff = FCcutoff)
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de/figures/kevin/Writeup1/Writeup1_de-volcano_wilcoxon.png"),
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
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de/figures/kevin/Writeup1/Writeup1_de-volcano_t.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")


