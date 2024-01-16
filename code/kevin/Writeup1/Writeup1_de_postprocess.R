rm(list=ls())

library(EnhancedVolcano)

load("../../../../out/kevin_demo/Demo1/Demo1_compute-de.RData")

dim(de_results)
head(de_results)
quantile(de_results$p_val_adj)

# determine p-value cutoff
tmp <- sort(de_results$p_val_adj, decreasing = F)
tmp <- tmp[tmp != 0]
pCutoff <- tmp[10]

# determine log fold change cutoff
FCcutoff <- sort(abs(de_results$avg_log2FC), decreasing = T)[100]

plot1 <-  EnhancedVolcano::EnhancedVolcano(de_results,
                                        lab = rownames(de_results),
                                        x = "avg_log2FC",
                                        y = "p_val_adj",
                                        pCutoff = pCutoff,
                                        FCcutoff = FCcutoff)
ggplot2::ggsave(filename = paste0("../../../../figures/kevin_demo/Demo1/Demo1_compute-de_volcano.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")
