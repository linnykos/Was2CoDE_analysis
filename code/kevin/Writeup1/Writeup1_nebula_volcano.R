rm(list=ls())
library(EnhancedVolcano)
set.seed(10)

load("~/lab/projects/subject-de/out/kevin/Writeup1/Writeup1_nebula.RData")

head(nebula_res$summary)
res <- nebula_res$summary

pval_vec <- res[,"p_CognitiveStatusNo dementia"]
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)
pCutoff <- max(pval_vec[idx])
FCcutoff <- quantile(abs(res[,"logFC_CognitiveStatusNo dementia"]), probs = 0.9)
xlim <- c(-1,1) * quantile(abs(res[,"logFC_CognitiveStatusNo dementia"]), probs = 0.99)

idx <- which(abs(res[,"logFC_CognitiveStatusNo dementia"]) <= max(xlim))
ylim <- c(0, max(-log10(pval_vec[idx])))

plot1 <-  EnhancedVolcano::EnhancedVolcano(
  res,
  lab = res$gene,
  x = "logFC_CognitiveStatusNo dementia",
  y = 'p_CognitiveStatusNo dementia',
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  xlim = xlim,
  ylim = ylim
)

ggplot2::ggsave(filename = "../../../figures/kevin/Writeup1/Writeup1_nebula_volcano.png",
                plot1, device = "png", width = 7, height = 7, units = "in")

