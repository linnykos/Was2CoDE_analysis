rm(list=ls())
library(EnhancedVolcano)
library(openxlsx)
library(ggplot2)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_eprs_nebula.RData")

num_levels <- length(nebula_res_list)
plot_list <- vector("list", length = num_levels)

for(lev in 1:num_levels){
  print(lev)
  res <- nebula_res_list[[lev]]$summary
  
  pval_vec <- res[,"p_ePRSLow"]
  pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
  idx <- which(pval_adj_vec <= 0.05)
  #define thresholds
  pCutoff <- max(pval_vec[idx])
  
  fc_vec <- res[,"logFC_ePRSLow"]
  fc_vec <- fc_vec[intersect(which(fc_vec < 6),
                             which(fc_vec > -6))]
  FCcutoff <- quantile(abs(fc_vec), probs = 0.9)
  xlim <- c(-1,1) * quantile(abs(fc_vec), probs = 0.99)
  
  idx <- which(abs(res[,"logFC_ePRSLow"]) <= max(xlim))
  ylim <- c(0, max(-log10(pval_vec[idx])))
  
  plot_list[[lev]] <- EnhancedVolcano::EnhancedVolcano(
    res,
    lab = res$gene,
    x = "logFC_ePRSLow",
    y = "p_ePRSLow",
    pCutoff = pCutoff,
    FCcutoff = FCcutoff,
    title = paste("NEBULA on cluster", lev-1),
    xlim = xlim,
    ylim = ylim
  )
}

grDevices::pdf(file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6/Writeup6_eprs_nebula_volcano.pdf"),
               width = 8, height = 8, onefile = TRUE)
for(lev in 1:num_levels){
  print(plot_list[[lev]])
}
grDevices::graphics.off()

#########################

for(lev in 1:num_levels){
  res <- nebula_res_list[[lev]]$summary
  
  df <- data.frame(
    p_val = res$p_ePRSLow,
    avg_log2FC = res$logFC_ePRSLow,
    p_val_adj = stats::p.adjust(res$p_ePRSLow, method = "BH")
  )
  rownames(df) <- res$gene
  df <- df[order(df$p_val_adj, decreasing = FALSE),]
  
  write.csv(df,
            file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/csv/kevin/Writeup6/Writeup6_eprs_nebula_cluster-", lev-1, ".csv"))
}
