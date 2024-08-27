rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(IdeasCustom)
library(ggplot2)
library(gridExtra)
library(ggrepel)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_Katie_Pseudobulk-DEseq2.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_Katie_was2_wilcox.RData")


gene_intersect_was2_deseq2 <- intersect(rownames(results_mat), 
                                        rownames(deseq2_res))

was2_pvalues <- results_mat[gene_intersect_was2_deseq2, "p_val"]
was2_logfc <- log2((results_mat[, "mean_dn"] ) / (results_mat[, "mean_nn"]+results_mat[, "mean_dd"])/2)[gene_intersect_was2_deseq2]
res <- data.frame(gene = gene_intersect_was2_deseq2,
                  was2_pvalue = was2_pvalues,
                  was2_logfc = was2_logfc,
                  deseq2_pvalue = deseq2_res[gene_intersect_was2_deseq2,"pvalue"],
                  deseq2_logfc = deseq2_res[gene_intersect_was2_deseq2,"log2FoldChange"])
rownames(res) <- res$gene

#define thresholds
pval_vec_deseq2 <- res[,"deseq2_pvalue"]
pval_adj_vec_deseq2 <- stats::p.adjust(pval_vec_deseq2, method = "BH")
idx_deseq2 <- which(pval_adj_vec_deseq2 <= 0.05)
pCutoff_deseq2 <- max(pval_vec_deseq2[idx_deseq2])

#define thresholds
pval_vec_was2 <- res[,"was2_pvalue"]
pval_adj_vec_was2 <- stats::p.adjust(pval_vec_was2, method = "BH")
idx_was2 <- which(pval_adj_vec_was2 <= 0.05)
pCutoff_was2 <- max(pval_vec_was2[idx_was2])

# Reorder the levels of SignificanceCategory
genes_above_threshold_deseq2 <- read.csv("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/genes_above_threshold_Pseudobulk-DEseq2.csv",
                                         header = FALSE, col.names = "gene")
genes_above_threshold_was2 <- res$gene[idx_was2]
significant_genes_was2 <- genes_above_threshold_was2
significant_genes_deseq2 <- genes_above_threshold_deseq2$gene

res$SignificanceCategory <- ifelse(rownames(res) %in% significant_genes_was2 & rownames(res) %in% significant_genes_deseq2, "Both",
                                   ifelse(rownames(res) %in% significant_genes_deseq2, "Pseudobulk",
                                          ifelse(rownames(res) %in% significant_genes_was2, "was2", "Other")))
res$SignificanceCategory <- factor(res$SignificanceCategory, levels = c("Other", "was2", "Pseudobulk","Both"))

res <- res[c(which(res$SignificanceCategory == "Other"),
             which(res$SignificanceCategory == "was2"),
             which(res$SignificanceCategory == "Pseudobulk"),
             which(res$SignificanceCategory == "Both")),
]
#print(res$gene[res$SignificanceCategory == "Both"])
# vector for significant genes
significant_genes <- res$SignificanceCategory != "Other"

############################
# Now create the plot to compare pvalues of was2 and Pseudobulk
correlation_pvalue_was2_deseq2 <- stats::cor(-log10(res$was2_pvalue),
                                               -log10(res$deseq2_pvalue))

plot0 <- ggplot(res, aes(x = -log10(was2_pvalue), y = -log10(deseq2_pvalue), color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "was2" = "red", "Pseudobulk" = "darkgreen", "Other" = "grey")) +
  geom_text(aes(label = ifelse(gene %in% c(genes_above_threshold_deseq2, genes_above_threshold_was2), as.character(gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  geom_text_repel(aes(label = ifelse(significant_genes, as.character(gene), "")),  # Conditionally label significant genes
                  box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = Inf) +
  labs(title = paste("Scatter Plot of PValues, Correlation =", round(correlation_pvalue_was2_deseq2, 2)), x = "-Log10 p-value (was2)", y = "-Log10 p-value (DESeq2)") +
  theme_minimal() +
  geom_vline(xintercept = -log10(pCutoff_was2), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -log10(pCutoff_deseq2), linetype = "dashed", color = "black")

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_PValues_was2-to-Pseudobulk-DEseq2.png",
                plot0, device = "png", width = 7, height = 7, units = "in")
############################

############################
# Compute thresholds for significant logFC based on the 90th percentile for both methods
FCcutoff_was2 <- quantile(abs(res$was2_logfc), probs = 0.9)
FCcutoff_deseq2 <- quantile(abs(res$deseq2_logfc), probs = 0.9)

# Set axis limits based on the larger 99th percentile of the absolute logFC values from both methods
xlims <- c(-1, 1) * max(quantile(abs(res$was2_logfc), probs = 0.99),
                        quantile(abs(res$deseq2_logfc), probs = 0.99))

# Compute correlation between logFC values, ensure no extra commas and properly handle NAs
correlation_logfc_was2_deseq2 <- stats::cor(res$was2_logfc, res$deseq2_logfc, use = "complete.obs")

# Plot
plot1 <- ggplot(res, aes(x = was2_logfc, y = deseq2_logfc, color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "was2" = "red", "Pseudobulk" = "darkgreen", "Other" = "grey")) +
  geom_text_repel(aes(label = ifelse(significant_genes, as.character(gene), "")),  # Conditionally label significant genes
                  box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = Inf) +
  labs(title = paste("Scatter Plot of Log2FC, Correlation =", round(correlation_logfc_was2_deseq2, 2)),
       x = "Log2 FC (was2)",
       y = "Log2 FC (deseq2)") +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + # Line of equality
  geom_vline(xintercept = c(-FCcutoff_was2, FCcutoff_was2), linetype = "dotted", color = "blue") +  
  geom_hline(yintercept = c(-FCcutoff_deseq2, FCcutoff_deseq2), linetype = "dotted", color = "purple") +  
  coord_cartesian(xlim = xlims, ylim = xlims)  # Set limits for x and y axes

ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Log2FC_Was2-decomp-to-Pseudobulk-DEseq2.png",
       plot1, device = "png", width = 7, height = 7, units = "in")

# Combine the plot0 and plot1
combined_plot_was2_to_deseq2 <- grid.arrange(plot0, plot1, ncol = 2)  # Arrange side by side
ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Was2-decomp-to-Pseudobulk-DEseq2.png",
       plot = combined_plot_was2_to_deseq2, device = "png", width = 14, height = 7, units = "in")
