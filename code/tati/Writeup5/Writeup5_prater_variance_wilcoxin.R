rm(list=ls())
library(ggplot2)
library(dplyr)
library(DESeq2)
library(Seurat)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(10)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

gene_vec <- Seurat::VariableFeatures(ss_data_norm[["RNA"]])
ss_data_norm <- subset(ss_data_norm, features = gene_vec)

extract_donor_class <- function(seurat_obj,
                                donor_variable_name,
                                class_variable_name,
                                class_variable_order){
  meta_data <- seurat_obj@meta.data
  stopifnot(all(c(class_variable_name, donor_variable_name) %in% colnames(meta_data)),
            length(class_variable_order) > 0,
            all(class_variable_order %in% meta_data[,class_variable_name]))
  
  tab_mat <- table(meta_data[,donor_variable_name], meta_data[,class_variable_name])
  condition <- apply(tab_mat, 1, function(x){colnames(tab_mat)[which(x!=0)]})
  condition <- factor(condition, levels = class_variable_order)
  
  return(condition)
}

condition <- extract_donor_class(
  seurat_obj = ss_data_norm,
  donor_variable_name = "Pt_ID",
  class_variable_name = "Study_Designation",
  class_variable_order = c("Ctrl", "AD")
)

# Variance
calculate_variance <- function(seurat_obj, donor_ids) {
  expression_matrix <- GetAssayData(seurat_obj, slot = "data")
  
  genes <- rownames(expression_matrix)
  donors <- unique(donor_ids)
  
  variance_matrix <- matrix(NA, nrow = length(donors), ncol = length(genes))
  colnames(variance_matrix) <- genes
  rownames(variance_matrix) <- donors
  
  for (gene in genes) {
    gene_expression <- expression_matrix[gene, ]
    for (donor in donors) {
      donor_cells <- gene_expression[donor_ids == donor]
      variance_matrix[as.character(donor), gene] <- var(donor_cells, na.rm = TRUE)
    }
  }
  
  return(variance_matrix)
}
meta_data <- ss_data_norm@meta.data
donor_ids <- meta_data$Pt_ID

variance_matrix <- calculate_variance(ss_data_norm, donor_ids)

#Wilcoxon test
perform_wilcoxon_test <- function(variance_matrix, condition) {
  stopifnot(nrow(variance_matrix) == length(condition), # check the number of donors
            all(sort(rownames(variance_matrix)) == sort(names(condition))), # check that the donor names match
            length(unique(names(condition))) == length(names(condition))) # check each donor has a unique name
  
  condition <- condition[rownames(variance_matrix)] # make the two objects in the same donor order
  
  genes <- colnames(variance_matrix)
  p_values <- sapply(genes, function(gene) {
    gene_variances <- variance_matrix[, gene]
    case_variances <- gene_variances[condition == levels(condition)[2]]
    control_variances <- gene_variances[condition == levels(condition)[1]]
    
    if (length(unique(case_variances)) > 1 & length(unique(control_variances)) > 1) {
      test_result <- wilcox.test(case_variances, control_variances)
      return(test_result$p.value)
    } else {
      return(NA)
    }
  })
  
  return(p_values)
}

wilcoxon_p_values <- perform_wilcoxon_test(variance_matrix, condition)

# some random things [[You can comment these lines out]]
######
quantile(wilcoxon_p_values)
sort(wilcoxon_p_values, decreasing = FALSE)[1:100]
wilcoxon_p_values_adjust <- stats::p.adjust(wilcoxon_p_values, method = "BH")
quantile(wilcoxon_p_values_adjust)
######

calculate_log2fc <- function(variance_matrix, condition) {
  stopifnot(nrow(variance_matrix) == length(condition), # check the number of donors
            all(sort(rownames(variance_matrix)) == sort(names(condition))), # check that the donor names match
            length(unique(names(condition))) == length(names(condition))) # check each donor has a unique name
  
  condition <- condition[rownames(variance_matrix)] # make the two objects in the same donor order
  
  genes <- colnames(variance_matrix)
  log2fc <- sapply(genes, function(gene) {
    gene_variances <- variance_matrix[, gene]
    case_variances <- gene_variances[condition == levels(condition)[2]]
    control_variances <- gene_variances[condition == levels(condition)[1]]
    
    mean_case <- mean(case_variances, na.rm = TRUE)
    mean_control <- mean(control_variances, na.rm = TRUE)
    
    if (is.na(mean_case) || is.na(mean_control) || mean_case == 0 || mean_control == 0) {
      return(NA)
    }
    
    return(log2(mean_case / mean_control))
  })
  
  return(log2fc)
}

log2fc <- calculate_log2fc(variance_matrix, condition)


result_df <- data.frame(Gene = colnames(variance_matrix), 
                        log2FC = log2fc, 
                        P_Value = wilcoxon_p_values)

# Remove NA values and adjust p-values
result_df <- result_df %>% 
  filter(!is.na(log2FC) & !is.na(P_Value)) %>%
  mutate(P_Value_adjusted = p.adjust(P_Value, method = "BH"),
         Significance = ifelse(P_Value_adjusted < 0.05, "Significant", "Not Significant"))

save(result_df, file = "~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_Prater_variance_results.RData")

# volcano plot
top_genes <- result_df %>%
  filter(Significance == "Significant") %>%
  top_n(10, wt = abs(log2FC))

volc_plot <- ggplot(result_df, aes(x = log2FC, y = -log10(P_Value), color = Significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3, box.padding = 0.5, point.padding = 0.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(title = "Volcano Plot of Gene Expression Variance Differences (Prater Data)",
       x = "Log2 Fold Change (Variance)",
       y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/Writeup5_Prater_variance_volcano.png",
       plot = volc_plot, device = "png", width = 7, height = 7, units = "in")

# Print significant genes
significant_genes <- result_df %>% 
  filter(Significance == "Significant") %>%
  arrange(P_Value)
print(significant_genes)

# GSEA analysis
teststat_vec <- result_df[,"log2FC"]
names(teststat_vec) <- result_df$Gene
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP",
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)

head(as.data.frame(gse))

gse_df <- as.data.frame(gse)
gse_df <- gse_df[gse_df$p.adjust <= 0.05,]
gse_df[1:100,"Description"]

# make ridgeplot of the gse results (https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html)

library(DOSE)
library(enrichplot)
library(ggplot2)

ridge_plot <- ridgeplot(gse) +
  labs(x = "Enrichment Score",
       y = "GO Biological Process",
       title = "Gene Set Enrichment Analysis Variance Prater's",
       subtitle = "Distribution of enrichment scores for top GO terms") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 8))

save_path <- "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/Writeup5_Prater_variance_GSEA_ridge_plot.png"
ggsave(save_path, ridge_plot, width = 12, height = 8, dpi = 300)


bar_plot <- barplot(gse, showCategory=20) +
  labs(title = "Top Enriched GO Terms",
       subtitle = "Gene Set Enrichment Analysis Results") +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(size = 8))


save_path <- "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/Writeup5_Prater_variance_GSEA_bar_plot.png"
ggsave(save_path, bar_plot, width = 12, height = 8, dpi = 300)

#dot_plot
plot1 <- enrichplot::dotplot(gse, showCategory=30) + ggplot2::ggtitle("dotplot for GSEA (ROSMAP, Variance test)")
plot1 <- plot1 + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/Writeup5_Prater_variance_GSEA_dot_plot.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

