library(ggplot2)
library(dplyr)
library(DESeq2)
library(Seurat)
library(ggrepel)
set.seed(10)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData")

# Extract metadata and expression data
meta_data <- seurat_obj@meta.data
seurat_obj <- GetAssayData(seurat_obj, layer = "data")  # Get normalized expression data
donor_labels <- meta_data$donor_id  # Donor ID for each cell
ADNC <- meta_data$ADNC  # Disease status: "Control" or "Case"
ADNC <- factor(ADNC, levels = c("Control", "Case"))

calculate_variance <- function(seurat_obj, donor_labels) {
  genes <- rownames(seurat_obj)
  donors <- unique(donor_labels)
  
  variance_matrix <- matrix(NA, nrow = length(donors), ncol = nrow(seurat_obj))
  colnames(variance_matrix) <- genes
  rownames(variance_matrix) <- donors
  
  for (gene_idx in 1:nrow(seurat_obj)) {
    for (donor in donors) {
      donor_cells <- seurat_obj[gene_idx, donor_labels == donor]
      variance_matrix[as.character(donor), gene_idx] <- var(donor_cells, na.rm = TRUE)
    }
  }
  
  return(variance_matrix)
}


variance_matrix <- calculate_variance(seurat_obj, donor_labels)

perform_wilcoxon_test <- function(variance_matrix, ADNC) {
  genes <- colnames(variance_matrix)
  p_values <- sapply(genes, function(gene) {
    gene_variances <- variance_matrix[, gene]
    case_variances <- gene_variances[ADNC == "Case"]
    control_variances <- gene_variances[ADNC == "Control"]
    
    if (length(unique(case_variances)) > 1 & length(unique(control_variances)) > 1) {
      test_result <- wilcox.test(case_variances, control_variances)
      return(test_result$p.value)
    } else {
      return(NA)
    }
  })
  
  return(p_values)
}

# Wilcoxon tests
wilcoxon_p_values <- perform_wilcoxon_test(variance_matrix, ADNC)

# initial dataframe
result_df <- data.frame(Gene = colnames(variance_matrix), P_Value = wilcoxon_p_values)
result_df <- result_df %>% filter(!is.na(P_Value))

calculate_log2fc <- function(variance_matrix, ADNC) {
  genes <- colnames(variance_matrix)
  log2fc <- sapply(genes, function(gene) {
    gene_variances <- variance_matrix[, gene]
    case_variances <- gene_variances[ADNC == "Case"]
    control_variances <- gene_variances[ADNC == "Control"]
    
    mean_case <- mean(case_variances, na.rm = TRUE)
    mean_control <- mean(control_variances, na.rm = TRUE)
    
    # Check for zero or NA means
    if (is.na(mean_case) || is.na(mean_control) || mean_case == 0 || mean_control == 0) {
      return(NA)
    }
    
    return(log2(mean_case / mean_control))
  })
  
  return(log2fc)
}


log2fc <- calculate_log2fc(variance_matrix, ADNC)
summary(log2fc)
sum(is.na(log2fc))

result_df$log2FC <- log2fc[result_df$Gene]
result_df <- result_df %>% filter(!is.na(log2FC)) # Remove rows with NA in log2FC
result_df$P_Value_adjusted <- stats::p.adjust(result_df$P_Value, method = "BH")
result_df$Significance <- ifelse(result_df$P_Value_adjusted < 0.05, "Significant", "Not Significant")

summary(result_df)

# Select genes
top_genes <- result_df %>%
  filter(Significance == "Significant") %>%
  top_n(10, wt = abs(log2FC))

volc_plot <- ggplot(result_df, aes(x = log2FC, y = -log10(P_Value), color = Significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3, box.padding = 0.5, point.padding = 0.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(title = "Volcano Plot of Gene Variance Differences",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup6/Writeup6_SEA-AD_variance_volcano.png",
       plot = volc_plot, device = "png", width = 7, height = 7, units = "in")

significant_genes <- result_df %>% 
  filter(Significance == "Significant") %>%
  arrange(P_Value)
print(significant_genes)

# Adjust plot limits?
# volc_plot <- volc_plot + 
#   xlim(-max(abs(result_df$log2FC)), max(abs(result_df$log2FC))) +
#   ylim(0, max(-log10(result_df$P_Value)))

#####################

# GSEA analysis
library(org.Hs.eg.db)
library(clusterProfiler)

teststat_vec <-result_df[,"log2FC"]
names(teststat_vec) <- rownames(result_df)
teststat_vec <- sort(teststat_vec, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

head(as.data.frame(gse))

