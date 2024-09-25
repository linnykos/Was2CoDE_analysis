rm(list=ls())
library(EnhancedVolcano)
library(openxlsx)
library(ggplot2)
library(IdeasCustom)
library(dplyr)
library(DESeq2)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_Pseudobulk-DEseq2.RData")

res <- deseq2_res

pval_vec <- res[,"pvalue"]
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)
#define thresholds
pCutoff <- max(pval_vec[idx])
FCcutoff <- quantile(abs(res[,"log2FoldChange"]), probs = 0.9)
xlim <- c(-1,1) * quantile(abs(res[,"log2FoldChange"]), probs = 0.99)

idx <- which(abs(res[,"log2FoldChange"]) <= max(xlim))
ylim <- c(0, max(-log10(pval_vec[idx])))

plot1 <- EnhancedVolcano::EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  xlim = xlim,
  ylim = ylim
)

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/Writeup7_ROSMAP_Pseudobulk-DEseq2_volcano.png",
                plot1, device = "png", width = 7, height = 7, units = "in")

#########################
res <- as.data.frame(deseq2_res)
res$gene <- rownames(res)
rownames(res) <- NULL

microglia_prater <- unique(microglia_prater_df$Gene[microglia_prater_df$DEG])
housekeeping_hounkpe <- housekeeping_hounkpe_df$Gene

#define thresholds
pval_vec <- res[,"pvalue"]
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)

pCutoff <- max(pval_vec[idx])
FCcutoff <- quantile(abs(res[,"log2FoldChange"]), probs = 0.9)
xlim <- c(-1,1) * quantile(abs(res[,"log2FoldChange"]), probs = 0.99)
ylim <- c(0, max(-log10(pval_vec[idx])))

res$GeneType <- ifelse(res$gene %in% microglia_prater, "Prater",
                       ifelse(res$gene %in% housekeeping_hounkpe, "Housekeeping", "Other"))


# First, reorder the levels of GeneType
res$GeneType <- factor(res$GeneType, levels = c("Other", "Housekeeping", "Prater"))
res <- res[c(which(res$GeneType == "Other"),
             which(res$GeneType == "Housekeeping"),
             which(res$GeneType == "Prater")),
]

# Now create the plot
plot0 <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = GeneType)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Prater" = "red", "Housekeeping" = "green", "Other" = "grey")) +
  geom_text(aes(label = ifelse(gene %in% c(microglia_prater, housekeeping_hounkpe), as.character(gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme_minimal() +
  xlim(xlim) + ylim(ylim) + 
  geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -log10(pCutoff), linetype = "dashed", color = "black")

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/Writeup7_ROSMAP_Pseudobulk-DEseq2_colored_by_GeneLists.png",
                plot0, device = "png", width = 7, height = 7, units = "in")

genes_above_threshold <- res %>%
  filter(-log10(pvalue) > -log10(pCutoff)) %>%
  select(gene)

print(genes_above_threshold)

write.table(genes_above_threshold, "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/Writeup7_ROSMAP_genes_above_threshold_Pseudobulk-DEseq2.csv", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
