rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(IdeasCustom)
library(ggplot2)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_Katie_Pseudobulk-DEseq2.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_nebula.RData")

set.seed(10)

gene_intersect <- intersect(rownames(deseq2_res), 
                            nebula_res$summary$gene)

res <- data.frame(gene = gene_intersect,
                  deseq2_pvalue = deseq2_res[gene_intersect,"pvalue"])
rownames(res) <- res$gene
tmp <- nebula_res$summary
rownames(tmp) <- tmp$gene
nebula_pvalue <- tmp[gene_intersect, "p_CognitiveStatusNo dementia"]
res$nebula_pvalue <- nebula_pvalue

microglia_prater <- unique(microglia_prater_df$Gene[microglia_prater_df$DEG])
housekeeping_hounkpe <- housekeeping_hounkpe_df$Gene

#define thresholds
pval_vec <- res[,"deseq2_pvalue"]
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)
pCutoff_deseq2 <- max(pval_vec[idx])

#define thresholds
pval_vec <- res[,"nebula_pvalue"]
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)
pCutoff_nebula <- max(pval_vec[idx])

res$GeneType <- ifelse(rownames(res) %in% microglia_prater, "Prater",
                       ifelse(rownames(res) %in% housekeeping_hounkpe, "Housekeeping", "Other"))

# First, reorder the levels of GeneType
res$GeneType <- factor(res$GeneType, levels = c("Other", "Housekeeping", "Prater"))
res <- res[c(which(res$GeneType == "Other"),
             which(res$GeneType == "Housekeeping"),
             which(res$GeneType == "Prater")),
]

# Now create the plot
correlation_value <- stats::cor(-log10(res$deseq2_pvalue),
                                -log10(res$nebula_pvalue))

plot0 <- ggplot(res, aes(x = -log10(deseq2_pvalue), y = -log10(nebula_pvalue), color = GeneType)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Prater" = "red", "Housekeeping" = "green", "Other" = "grey")) +
  geom_text(aes(label = ifelse(gene %in% c(microglia_prater, housekeeping_hounkpe), as.character(gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  labs(title = paste("Volcano Plot, Correlation =", round(correlation_value, 2)), x = "-Log10 p-value (DESeq2)", y = "-Log10 p-value (Nebula)") +
  theme_minimal() +
  geom_vline(xintercept = -log10(pCutoff_deseq2), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -log10(pCutoff_nebula), linetype = "dashed", color = "black")

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Nebula-to-Pseudobulk-DEseq2.png",
                plot0, device = "png", width = 7, height = 7, units = "in")
