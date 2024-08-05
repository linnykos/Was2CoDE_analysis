rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(was2)
library(IdeasCustom)
library(ggplot2)
library(gridExtra)
library(ggrepel)
set.seed(10)
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_prater_esvd.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_Katie_was2_wilcox.RData")

# Intersect genes from esvd and WAS2 datasets
gene_intersect_esvd_was2 <- intersect(names(eSVD_obj$teststat_vec), rownames(results_mat))
was2_pvalues <- results_mat[gene_intersect_esvd_was2, "p_val"]
was2_logfc <- log2((results_mat[, "mean_dn"] ) / (results_mat[, "mean_nn"]+results_mat[, "mean_dd"])/2)[gene_intersect_esvd_was2]

esvd_p_values <- 10^(-eSVD_obj$pvalue_list$log10pvalue[gene_intersect_esvd_was2])
esvd_logfc <- (log2(eSVD_obj$case_mean) - log2(eSVD_obj$control))[gene_intersect_esvd_was2]


res <- data.frame(gene = gene_intersect_esvd_was2,
                  was2_pvalue = was2_pvalues,
                  was2_logfc = was2_logfc,
                  esvd_pvalue = esvd_p_values,
                  esvd_logfc = esvd_logfc)

#define thresholds
pval_vec_was2 <- res[,"was2_pvalue"]
pval_adj_vec_was2 <- stats::p.adjust(pval_vec_was2, method = "BH")
idx_was2 <- which(pval_adj_vec_was2 <= 0.05)
pCutoff_was2 <- max(pval_vec_was2[idx_was2])

#define thresholds
pval_vec_esvd <- res[,"esvd_pvalue"]
pval_adj_vec_esvd <- stats::p.adjust(pval_vec_esvd, method = "BH")
idx_esvd <- which(pval_adj_vec_esvd <= 0.05)
pCutoff_esvd <- max(pval_vec_esvd[idx_esvd])

############################
# # Reorder the levels of GeneType
############################
# microglia_prater <- unique(microglia_prater_df$Gene[microglia_prater_df$DEG])
# housekeeping_hounkpe <- housekeeping_hounkpe_df$Gene
# 
# res$GeneType <- ifelse(rownames(res) %in% microglia_prater, "Prater",
#                        ifelse(rownames(res) %in% housekeeping_hounkpe, "Housekeeping", "Other"))
# 
# res$GeneType <- factor(res$GeneType, levels = c("Other", "Housekeeping", "Prater"))
# 
# res <- res[c(which(res$GeneType == "Other"),
#              which(res$GeneType == "Housekeeping"),
#              which(res$GeneType == "Prater")),
# ]
############################
############################
# Reorder the levels of SignificanceCategory
genes_above_threshold_eSVD <- read.csv("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/genes_above_threshold_eSVD.csv",
                                       header = FALSE, col.names = "gene")

# Extract gene names based on significant indices
genes_above_threshold_was2 <- res$gene[idx_was2]
# Create vectors from the data frames for easier checking
significant_genes_was2 <- genes_above_threshold_was2
significant_genes_esvd <- genes_above_threshold_eSVD$gene

res$SignificanceCategory <- ifelse(rownames(res) %in% significant_genes_was2 & rownames(res) %in% significant_genes_esvd, "Both",
                                   ifelse(rownames(res) %in% significant_genes_was2, "was2",
                                          ifelse(rownames(res) %in% significant_genes_esvd, "esvd", "Other")))
res$SignificanceCategory <- factor(res$SignificanceCategory, levels = c("Other", "esvd", "was2","Both"))

res <- res[c(which(res$SignificanceCategory == "Other"),
             which(res$SignificanceCategory == "esvd"),
             which(res$SignificanceCategory == "was2"),
             which(res$SignificanceCategory == "Both")),
]
# vector for significant genes
significant_genes <- res$SignificanceCategory != "Other"

############################


############################

# Create the plot to compare pvalues of was2 and esvd
correlation_pvalue_was2_esvd <- stats::cor(-log10(res$was2_pvalue),
                                             -log10(res$esvd_pvalue))

plot0 <- ggplot(res, aes(x = -log10(was2_pvalue), y = -log10(esvd_pvalue), color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "was2" = "red", "esvd" = "darkgreen", "Other" = "grey")) +
  geom_text(aes(label = ifelse(gene %in% c(genes_above_threshold_was2, genes_above_threshold_eSVD), as.character(gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  geom_text_repel(aes(label = ifelse(significant_genes, as.character(gene), "")),  # Conditionally label significant genes
                  box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = Inf) +
  labs(title = paste("Scatter Plot of PValues, Correlation =", round(correlation_pvalue_was2_esvd, 2)), x = "-Log10 p-value (was2)", y = "-Log10 p-value (esvd)") +
  theme_minimal() +
  geom_vline(xintercept = -log10(pCutoff_was2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pCutoff_esvd), linetype = "dashed", color = "black")

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_PValues_Was2-decomp-to-eSVD.png",
                plot0, device = "png", width = 7, height = 7, units = "in")
############################

############################
# Create the plot to compare logfc of was2 and esvd

# Compute thresholds for significant logFC based on the 90th percentile for both methods
FCcutoff_was2 <- quantile(abs(res$was2_logfc), probs = 0.9)
FCcutoff_esvd <- quantile(abs(res$esvd_logfc), probs = 0.9)

# Set axis limits based on the larger 99th percentile of the absolute logFC values from both methods
xlims <- c(-1, 1) * max(quantile(abs(res$was2_logfc), probs = 0.99),
                        quantile(abs(res$esvd_logfc), probs = 0.99))

# Compute correlation between logFC values
correlation_logfc_was2_esvd <- stats::cor(res$was2_logfc, res$esvd_logfc)

# Plot
plot1 <- ggplot(res, aes(x = was2_logfc, y = esvd_logfc, color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "was2" = "red", "esvd" = "darkgreen", "Other" = "grey")) +
  labs(title = paste("Scatter Plot of Log2FC, Correlation =", round(correlation_logfc_was2_esvd, 2)),
       x = "Log2 FC (was2)",
       y = "Log2 FC (esvd)") +
  geom_text_repel(aes(label = ifelse(significant_genes, as.character(gene), "")),  # Conditionally label significant genes
                  box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = Inf) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + # Line of equality
  geom_vline(xintercept = c(-FCcutoff_was2, FCcutoff_was2), linetype = "dotted", color = "blue") +  # Vertical lines for was2 FC cutoff
  geom_hline(yintercept = c(-FCcutoff_esvd, FCcutoff_esvd), linetype = "dotted", color = "purple") +  # Horizontal lines for esvd FC cutoff
  coord_cartesian(xlim = xlims, ylim = xlims)  # Set limits for x and y axes

ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Log2FC_Was2-decomp-to-eSVD.png",
       plot1, device = "png", width = 7, height = 7, units = "in")

# Combine the plot0 and plot1
combined_plot_was2_to_esvd <- grid.arrange(plot0, plot1, ncol = 2)  # Arrange side by side
ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Was2-decomp-to-eSVD.png",
       plot = combined_plot_was2_to_esvd, device = "png", width = 14, height = 7, units = "in")
