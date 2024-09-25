rm(list=ls())
library(EnhancedVolcano)
library(openxlsx)
library(ggplot2)
library(IdeasCustom)
library(dplyr)
library(DESeq2)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_eSVD.RData")
# Extract components from eSVD_obj
gene_names <- names(eSVD_obj$teststat_vec)
log2FoldChange <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
log10_pvalues <- eSVD_obj$pvalue_list$log10pvalue
p_values <- 10^(-log10_pvalues)
pval_adj_vec <- p.adjust(p_values, method = "BH")

# Define thresholds
significant_threshold <- 0.05
idx <- which(pval_adj_vec <= significant_threshold)
pCutoff <- max(p_values[idx])
FCcutoff <- quantile(abs(log2FoldChange), probs = 0.9)
xlim <- c(-1,1) * quantile(abs(log2FoldChange), probs = 0.99)
ylim <- c(0, max(-log10(p_values[idx])))

# Create a data frame for the volcano plot
res <- data.frame(Gene = gene_names, log2FoldChange = log2FoldChange, PValue = p_values)

microglia_prater <- unique(microglia_prater_df$Gene[microglia_prater_df$DEG])
housekeeping_hounkpe <- housekeeping_hounkpe_df$Gene
res$GeneType <- ifelse(res$Gene %in% microglia_prater, "Prater",
                       ifelse(res$Gene %in% housekeeping_hounkpe, "Housekeeping", "Other"))

# Reorder the levels of GeneType
res$GeneType <- factor(res$GeneType, levels = c("Other", "Housekeeping", "Prater"))
res <- res[c(which(res$GeneType == "Other"),
             which(res$GeneType == "Housekeeping"),
             which(res$GeneType == "Prater")),
]

# Create the volcano plot
plot0 <- ggplot(res, aes(x = log2FoldChange, y = -log10(PValue), color = GeneType)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Prater" = "red", "Housekeeping" = "green", "Other" = "grey")) +
  geom_text(aes(label = ifelse(Gene %in% c(microglia_prater, housekeeping_hounkpe), as.character(Gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme_minimal() +
  xlim(xlim) + ylim(ylim) + 
  geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -log10(pCutoff), linetype = "dashed", color = "black")

ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/Writeup7_ROSMAP_eSVD_volcano_colored_by_GeneLists.png",
       plot0, device = "png", width = 7, height = 7, units = "in")

# genes above significance threshold
genes_above_threshold <- res[res$PValue < pCutoff, ]
print(genes_above_threshold$Gene)
write.table(genes_above_threshold$Gene, "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/genes_above_threshold_eSVD.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
