rm(list=ls())
library(EnhancedVolcano)
set.seed(10)

load("~/lab/projects/subject-de/out/kevin/Writeup1/Writeup1_nebula.RData")

head(nebula_res$summary)
res <- nebula_res$summary

plot1 <-  EnhancedVolcano::EnhancedVolcano(res,
                                           lab = rownames(res),
                                           x = 'log2FoldChange',
                                           y = 'pvalue')