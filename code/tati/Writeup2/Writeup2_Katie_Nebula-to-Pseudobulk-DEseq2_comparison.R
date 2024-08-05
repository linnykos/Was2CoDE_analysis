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
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_nebula.RData")


gene_intersect_deseq2_nebula <- intersect(rownames(deseq2_res), 
                                          nebula_res$summary$gene)

res <- data.frame(gene = gene_intersect_deseq2_nebula,
                  deseq2_pvalue = deseq2_res[gene_intersect_deseq2_nebula,"pvalue"],
                  deseq2_logfc = deseq2_res[gene_intersect_deseq2_nebula,"log2FoldChange"])
rownames(res) <- res$gene
tmp <- nebula_res$summary
rownames(tmp) <- tmp$gene
nebula_pvalue <- tmp[gene_intersect_deseq2_nebula, "p_CognitiveStatusNo dementia"]
res$nebula_pvalue <- nebula_pvalue
nebula_logfc <- tmp[gene_intersect_deseq2_nebula, "logFC_CognitiveStatusNo dementia"]
res$nebula_logfc <-nebula_logfc 

# microglia_prater <- unique(microglia_prater_df$Gene[microglia_prater_df$DEG])
# housekeeping_hounkpe <- housekeeping_hounkpe_df$Gene
# res$GeneType <- ifelse(rownames(res) %in% microglia_prater, "Prater",
#                        ifelse(rownames(res) %in% housekeeping_hounkpe, "Housekeeping", "Other"))
#reorder the levels of GeneType
# res$GeneType <- factor(res$GeneType, levels = c("Other", "Housekeeping", "Prater"))
# res <- res[c(which(res$GeneType == "Other"),
#              which(res$GeneType == "Housekeeping"),
#              which(res$GeneType == "Prater")),
# ]

#define thresholds
pval_vec_deseq2 <- res[,"deseq2_pvalue"]
pval_adj_vec_deseq2 <- stats::p.adjust(pval_vec_deseq2, method = "BH")
idx_deseq2 <- which(pval_adj_vec_deseq2 <= 0.05)
pCutoff_deseq2 <- max(pval_vec_deseq2[idx_deseq2])

#define thresholds
pval_vec_nebula <- res[,"nebula_pvalue"]
pval_adj_vec_nebula <- stats::p.adjust(pval_vec_nebula, method = "BH")
idx_nebula <- which(pval_adj_vec_nebula <= 0.05)
pCutoff_nebula <- max(pval_vec_nebula[idx_nebula])

# Reorder the levels of SignificanceCategory
genes_above_threshold_deseq2 <- read.csv("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/genes_above_threshold_Pseudobulk-DEseq2.csv",
                                       header = FALSE, col.names = "gene")
genes_above_threshold_nebula <- read.csv("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/genes_above_threshold_nebula.csv",
                                         header = FALSE, col.names = "gene")
# Create vectors from the data frames for easier checking
significant_genes_deseq2 <- genes_above_threshold_deseq2$gene
significant_genes_nebula <- genes_above_threshold_nebula$gene

res$SignificanceCategory <- ifelse(rownames(res) %in% significant_genes_deseq2 & rownames(res) %in% significant_genes_nebula, "Both",
                                   ifelse(rownames(res) %in% significant_genes_deseq2, "Pseudobulk",
                                          ifelse(rownames(res) %in% significant_genes_nebula, "NEBULA", "Other")))
res$SignificanceCategory <- factor(res$SignificanceCategory, levels = c("Other", "NEBULA", "Pseudobulk","Both"))

res <- res[c(which(res$SignificanceCategory == "Other"),
             which(res$SignificanceCategory == "NEBULA"),
             which(res$SignificanceCategory == "Pseudobulk"),
             which(res$SignificanceCategory == "Both")),
]

# vector for significant genes
significant_genes <- res$SignificanceCategory != "Other"
############################


############################
# Now create the plot to compare pvalues of Pseudobulk and NEBULA
correlation_pvalue_deseq2_nebula <- stats::cor(-log10(res$deseq2_pvalue),
                                -log10(res$nebula_pvalue))

plot0 <- ggplot(res, aes(x = -log10(deseq2_pvalue), y = -log10(nebula_pvalue), color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "Pseudobulk" = "red", "NEBULA" = "darkgreen", "Other" = "grey")) +
  geom_text_repel(aes(label = ifelse(significant_genes, as.character(gene), "")),  # Conditionally label significant genes
                  box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = Inf) +
  labs(title = paste("Scatter Plot of PValues, Correlation =", round(correlation_pvalue_deseq2_nebula, 2)), x = "-Log10 p-value (DESeq2)", y = "-Log10 p-value (Nebula)") +
  theme_minimal() +
  geom_vline(xintercept = -log10(pCutoff_deseq2), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -log10(pCutoff_nebula), linetype = "dashed", color = "black")


ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_PValues_Nebula-to-Pseudobulk-DEseq2.png",
                plot0, device = "png", width = 7, height = 7, units = "in")
############################

############################
# Compute thresholds for significant logFC based on the 90th percentile for both methods
FCcutoff_deseq2 <- quantile(abs(res$deseq2_logfc), probs = 0.9)
FCcutoff_nebula <- quantile(abs(res$nebula_logfc), probs = 0.9)

# Set axis limits based on the larger 99th percentile of the absolute logFC values from both methods
xlims <- c(-1, 1) * max(quantile(abs(res$deseq2_logfc), probs = 0.99),
                        quantile(abs(res$nebula_logfc), probs = 0.99))

# Compute correlation between logFC values
correlation_logfc_deseq2_nebula <- stats::cor(res$deseq2_logfc, res$nebula_logfc)

# Plot
plot1 <- ggplot(res, aes(x = deseq2_logfc, y = nebula_logfc, color = SignificanceCategory)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Both" = "purple", "Pseudobulk" = "red", "NEBULA" = "darkgreen", "Other" = "grey")) +
  geom_text_repel(aes(label = ifelse(significant_genes, as.character(gene), "")),  # Conditionally label significant genes
                  box.padding = 0.5, point.padding = 0.3, size = 3, max.overlaps = Inf) +
  labs(title = paste("Scatter Plot of Log2FC, Correlation =", round(correlation_logfc_deseq2_nebula, 2)),
       x = "Log2 FC (deseq2)",
       y = "Log2 FC (Nebula)") +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + # Line of equality
  geom_vline(xintercept = c(-FCcutoff_deseq2, FCcutoff_deseq2), linetype = "dotted", color = "blue") +  # Vertical lines for deseq2 FC cutoff
  geom_hline(yintercept = c(-FCcutoff_nebula, FCcutoff_nebula), linetype = "dotted", color = "purple") +  # Horizontal lines for Nebula FC cutoff
  coord_cartesian(xlim = xlims, ylim = xlims)  # Set limits for x and y axes

ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Log2FC_Pseudobulk-DEseq2-to-NEBULA.png",
       plot1, device = "png", width = 7, height = 7, units = "in")

# Combine the plot0 and plot1
combined_plot_deseq2_to_nebula <- grid.arrange(plot0, plot1, ncol = 2)  # Arrange side by side
ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_Pseudobulk-DEseq2-to-NEBULA.png",
       plot = combined_plot_deseq2_to_nebula, device = "png", width = 14, height = 7, units = "in")
