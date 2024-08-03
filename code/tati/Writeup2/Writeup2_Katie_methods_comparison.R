rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(IdeasCustom)
library(ggplot2)
library(gridExtra)
library(ggrepel)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_nebula.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_prater_esvd.RData")
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_Katie_Pseudobulk-DEseq2.RData")

gene_intersect_esvd_nebula <- intersect(names(eSVD_obj$teststat_vec), 
                            nebula_res$summary$gene)
esvd_p_values <- 10^(-eSVD_obj$pvalue_list$log10pvalue[gene_intersect_esvd_nebula])
esvd_logfc <- (log2(eSVD_obj$case_mean) - log2(eSVD_obj$control))[gene_intersect_esvd_nebula]
res <- data.frame(gene = gene_intersect_esvd_nebula,
                  esvd_pvalue = esvd_p_values,
                  esvd_logfc = esvd_logfc)
rownames(res) <- res$gene
tmp <- nebula_res$summary
rownames(tmp) <- tmp$gene
nebula_pvalue <- tmp[gene_intersect_esvd_nebula, "p_CognitiveStatusNo dementia"]
res$nebula_pvalue <- nebula_pvalue
nebula_logfc <- tmp[gene_intersect_esvd_nebula, "logFC_CognitiveStatusNo dementia"]
res$nebula_logfc <-nebula_logfc 

#define thresholds
pval_vec_esvd <- res[,"esvd_pvalue"]
pval_adj_vec_esvd <- stats::p.adjust(pval_vec_esvd, method = "BH")
idx_esvd <- which(pval_adj_vec_esvd <= 0.05)
pCutoff_esvd <- max(pval_vec_esvd[idx_esvd])

#define thresholds
pval_vec_nebula <- res[,"nebula_pvalue"]
pval_adj_vec_nebula <- stats::p.adjust(pval_vec_nebula, method = "BH")
idx_nebula <- which(pval_adj_vec_nebula <= 0.05)
pCutoff_nebula <- max(pval_vec_nebula[idx_nebula])
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
genes_above_threshold_nebula <- read.csv("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/genes_above_threshold_nebula.csv",
                                         header = FALSE, col.names = "gene")
# Create vectors from the data frames for easier checking
significant_genes_eSVD <- genes_above_threshold_eSVD$gene
significant_genes_nebula <- genes_above_threshold_nebula$gene

res$SignificanceCategory <- ifelse(rownames(res) %in% significant_genes_eSVD & rownames(res) %in% significant_genes_nebula, "Both",
                                   ifelse(rownames(res) %in% significant_genes_eSVD, "eSVD",
                                          ifelse(rownames(res) %in% significant_genes_nebula, "NEBULA", "Other")))
res$SignificanceCategory <- factor(res$SignificanceCategory, levels = c("Other", "NEBULA", "eSVD","Both"))

res <- res[c(which(res$SignificanceCategory == "Other"),
             which(res$SignificanceCategory == "NEBULA"),
             which(res$SignificanceCategory == "eSVD"),
             which(res$SignificanceCategory == "Both")),
]
############################


############################

# Create the plot to compare pvalues of esvd and nebula
correlation_pvalue_esvd_nebula <- stats::cor(-log10(res$esvd_pvalue),
                                             -log10(res$nebula_pvalue))

plot0 <- ggplot(res, aes(x = -log10(esvd_pvalue), y = -log10(nebula_pvalue), color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "eSVD" = "red", "NEBULA" = "darkgreen", "Other" = "grey")) +
  geom_text(aes(label = ifelse(gene %in% c(genes_above_threshold_eSVD, genes_above_threshold_nebula), as.character(gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  labs(title = paste("Scatter Plot of PValues, Correlation =", round(correlation_pvalue_esvd_nebula, 2)), x = "-Log10 p-value (eSVD)", y = "-Log10 p-value (NEBULA)") +
  theme_minimal() +
  geom_vline(xintercept = -log10(pCutoff_esvd), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pCutoff_nebula), linetype = "dashed", color = "black")

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_PValues_eSVD-to-NEBULA.png",
                plot0, device = "png", width = 7, height = 7, units = "in")
############################

############################
# Create the plot to compare logfc of esvd and nebula

# Compute thresholds for significant logFC based on the 90th percentile for both methods
FCcutoff_esvd <- quantile(abs(res$esvd_logfc), probs = 0.9)
FCcutoff_nebula <- quantile(abs(res$nebula_logfc), probs = 0.9)

# Set axis limits based on the larger 99th percentile of the absolute logFC values from both methods
xlims <- c(-1, 1) * max(quantile(abs(res$esvd_logfc), probs = 0.99),
                        quantile(abs(res$nebula_logfc), probs = 0.99))

# Compute correlation between logFC values
correlation_logfc_esvd_nebula <- stats::cor(res$esvd_logfc, res$nebula_logfc)

# Plot
plot1 <- ggplot(res, aes(x = esvd_logfc, y = nebula_logfc, color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "eSVD" = "red", "NEBULA" = "darkgreen", "Other" = "grey")) +
  geom_text(aes(label = ifelse(gene %in% c(microglia_prater, housekeeping_hounkpe), as.character(gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  labs(title = paste("Scatter Plot of Log2FC, Correlation =", round(correlation_logfc_esvd_nebula, 2)),
       x = "Log2 FC (eSVD)",
       y = "Log2 FC (Nebula)") +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + # Line of equality
  geom_vline(xintercept = c(-FCcutoff_esvd, FCcutoff_esvd), linetype = "dotted", color = "blue") +  # Vertical lines for eSVD FC cutoff
  geom_hline(yintercept = c(-FCcutoff_nebula, FCcutoff_nebula), linetype = "dotted", color = "purple") +  # Horizontal lines for Nebula FC cutoff
  coord_cartesian(xlim = xlims, ylim = xlims)  # Set limits for x and y axes

ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Log2FC_eSVD-to-NEBULA.png",
       plot1, device = "png", width = 7, height = 7, units = "in")

# Combine the plot0 and plot1
combined_plot_esvd_to_nebula <- grid.arrange(plot0, plot1, ncol = 2)  # Arrange side by side
ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_eSVD-to-NEBULA.png",
       plot = combined_plot_esvd_to_nebula, device = "png", width = 14, height = 7, units = "in")
