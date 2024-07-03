# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData")
data_matrix <- SeuratObject::LayerData(ss_data_norm, layer = "data", assay = "RNA") # genes by cells
# dense_matrix <- as.matrix(data_matrix)
# genes <- rownames(dense_matrix)
# indices <- colnames(dense_matrix)
# dim(dense_matrix)

# first, create a list (one element per donor) that contains the indicies of all the cells for that donor

# then, assign colors to each donor (red = AD, blue = control)

#####

# iterate over a list of genes: (Grab this from your eSVD-DE or NEBULA analyses)
# these give you specific rows of your data matrix to work with

## for each gene, do the following:

## first, find the correct row for that gene in the data matrix
## then, create a data frame (rows = cells, and the columns can be: 
### 1) the gene expression value,
### 2) the donor that cell comes from
### 3) the status of the donor (AD or not)
## pass that data frame into ggplot to make the boxplots

meta_data <- ss_data_norm@meta.data
patient_ids <- meta_data$Pt_ID

# Create a data frame to map indices to patient IDs
df <- data.frame(Indices = indices, Pt_ID = patient_ids[match(indices, rownames(meta_data))])
head(df)

expression_data <- as.data.frame(dense_matrix) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "Indices", values_to = "Expression")
head(expression_data)
dim(expression_data)

# Merge expression data with df to get Pt_ID
expression_data <- expression_data %>%
  left_join(df, by = "Indices")


genes_of_interest <- c("SORCS1", "DPYS")

gene_indices <- which(rownames(data) %in% genes_of_interest)
gene_expression <- data[gene_indices, ]
gene_expression_df <- as.data.frame(as.matrix(gene_expression))
colnames(gene_expression_df) <- genes_of_interest

gene_expression <- gene_expression %>%
  rownames_to_column(var = "Cell") %>%
  gather(key = "Gene", value = "Expression", -Cell)


gene_expression <- gene_expression %>%
  mutate(Patient_ID = patient_ids[match(Cell, cell_names)])
# Merge the expression data with the patient IDs
gene_expression_id <- gene_expression %>%
  left_join(meta_data, by = "Cell")

for (gene in genes_of_interest) {
  expression_data <- FetchData(ss_data_norm, vars = c("Pt_ID", gene))
  plot <- ggplot(expression_data, aes(x = Pt_ID, y = !!sym(gene), fill = Pt_ID)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    theme_ipsum() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ) +
    ggtitle(paste("Boxplot for Gene:", gene)) +
    xlab("Donors") +
    ylab("Expression")
  
  # Save the plot
  ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_", gene, "_boxplot.png"),
         plot = plot, device = "png", width = 5, height = 7, units = "in")
}
  #########################################
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_nebula.RData")
nebula_summary <- nebula_res$summary

significant_genes <- c("SORCS1", "DPYS", "ASTN1", "GLT1D1", "BRINP2", "PBX1", "MAPK10", "RP11-762E8.1", "FSTL5", "ARHGAP10")

filtered_data <- nebula_summary %>%
  filter(gene %in% significant_genes)

for (gene in significant_genes) {
  # Filter the data for the current gene
  gene_data <- nebula_summary %>%
    filter(gene == !!gene)
  
  plot1 <- ggplot(gene_data, aes(x = `logFC_CognitiveStatusNo dementia`, y = `logFC_(Intercept)`, fill = `logFC_CognitiveStatusNo dementia`)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    theme_ipsum() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ) +
    ggtitle(paste("Boxplot for Gene:", gene)) +
    xlab("Cognitive Status (AD vs. non-AD)") +
    ylab("logFC (Intercept)")

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_nebula_boxplot_", gene, ".png"),
                plot1, device = "png", width = 5, height = 7, units = "in")
}