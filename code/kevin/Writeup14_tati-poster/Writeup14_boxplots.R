# from https://github.com/linnykos/subject-de/blob/tati/code/tati/Writeup5/Writeup5_boxplots.R
rm(list=ls())
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)
library(ggrastr)
library(ggplot2)

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup14/"

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

count_mat <- SeuratObject::LayerData(ss_data_norm, 
                                     layer = "data", 
                                     assay = "RNA")

# first, create a list (one element per donor) that contains the indicies of all the cells for that donor
meta_data <- ss_data_norm@meta.data
Pt_ids <- meta_data$Pt_ID
cell_names <- colnames(count_mat)
cell_df <- data.frame(cell_name = cell_names, Pt_ID = Pt_ids, stringsAsFactors = FALSE)

color_mapping <- meta_data %>%
  select(Pt_ID, Study_Designation) %>%
  distinct() %>%
  mutate(color = ifelse(Study_Designation == "AD", "#6a0dad", "#f4c20d"))

# named vector
color_mapping <- setNames(color_mapping$color, color_mapping$Pt_ID)

# Find genes with adjusted p-values <= 0.05
# genes_of_interest <- names(pvalue_vec_adjusted)[pvalue_vec_adjusted <= 0.05]

genes_of_interest <- c("ANXA2")


pdf(file = paste0(plot_folder, "Writeup14_boxplots.pdf"),
    width = 6, height = 6, onefile = TRUE)

# Iterate over each gene of interest and create a boxplot
for (gene in genes_of_interest) {
  # Find the correct row for that gene in the count matrix
  gene_index <- which(rownames(count_mat) == gene)
  
  # Extract expression for that gene
  expression_value <- count_mat[gene_index, ]
  
  df <- data.frame(
    expression = expression_value,
    Pt_ID = meta_data$Pt_ID,
    Study_Designation = meta_data$Study_Designation
  )
  
  # ylim <- c(min(df$expression), stats::quantile(df$expression, prob = 0.95))
  
  # Set ordering of donors
  ordering <- c(names(color_mapping)[color_mapping == "#6a0dad"],
                names(color_mapping)[color_mapping == "#f4c20d"])
  df$Pt_ID <- factor(df$Pt_ID, levels = ordering)
  
  # Generate boxplot with was2 LFC labeling
  plot <- ggplot(df, aes(x = Pt_ID, y = expression, fill = Pt_ID)) +
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = color_mapping) +
    theme_minimal() +
    ylim(c(0, 0.05)) +
    labs(title = paste("Expression of", gene), 
         x = "Donor", 
         y = "Expression") +
    theme(axis.text.x = element_blank()) +
    Seurat::NoLegend()
  
  # Print the plot to the PDF file
  print(plot)
}

# Close the PDF device
graphics.off()