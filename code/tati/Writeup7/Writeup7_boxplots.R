# Libraries
rm(list=ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup11/Writeup11_rosmap_scVI-postprocessed.RData") 

# Extract the count matrix (genes by cells)
count_mat <- SeuratObject::LayerData(seurat_obj, layer = "data", assay = "RNA")

# Get the patient IDs and create a data frame for cells and their corresponding donor IDs
Pt_ID <- seurat_obj$Pt_ID
cell_names <- colnames(count_mat)
cell_df <- data.frame(cell_name = cell_names, Pt_ID = Pt_ID, stringsAsFactors = FALSE)

# Extract the metadata from the Seurat object as a data frame
metadata_df <- seurat_obj@meta.data

# Assign colors to each donor (red = AD, blue = control)
color_mapping <- metadata_df %>%
  select(Pt_ID, ADpath) %>%
  distinct() %>%
  mutate(color = ifelse(ADpath == "yes", "red", "blue"))

# Inspect the color_mapping data frame
print(color_mapping)

# named vector
color_mapping <- setNames(color_mapping$color, color_mapping$Pt_ID)

#####
load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_eSVD.RData")
pvalue_vec <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
names(pvalue_vec) <- names(eSVD_obj$teststat_vec)
pvalue_vec_adjusted <- stats::p.adjust(pvalue_vec, method = "BH")
genes_of_interest <- names(pvalue_vec_adjusted)[pvalue_vec_adjusted <= 0.05]
genes_of_interest <- intersect(genes_of_interest, rownames(count_mat))
genes_of_interest

# Create PDF for boxplots
grDevices::pdf(file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/Writeup7_eSVD_boxplots.pdf"),
               width = 8, height = 8, onefile = TRUE)
# Iterate over each gene and generate the boxplot
for (gene in genes_of_interest) {
  tryCatch({
    # Find the correct row for that gene in the count matrix
    gene_index <- which(rownames(count_mat) == gene)
    
    # Extract expression for that gene across all cells
    expression_value <- count_mat[gene_index, ]
    
    # Create a data frame with expression values and metadata (Pt_ID, ADpath)
    df <- data.frame(
      expression = expression_value,
      Pt_ID = metadata_df$Pt_ID,  # Using the metadata from the Seurat object
      ADpath = metadata_df$ADpath  # AD status
    )
    
    # Apply the color mapping based on the AD status (red = AD, blue = control)
    df$color <- color_mapping[df$Pt_ID]
    
    # Order the patients by their AD status (red first, then blue)
    ordering <- unique(df$Pt_ID[order(df$color == "red", decreasing = TRUE)])
    df$Pt_ID <- factor(df$Pt_ID, levels = ordering)
    
    # Create the boxplot
    plot <- ggplot(df, aes(x = Pt_ID, y = expression, fill = Pt_ID)) +
      geom_boxplot(outlier.shape = NA) +  # Avoid showing outliers
      scale_fill_manual(values = setNames(df$color, df$Pt_ID)) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = paste("Expression of", gene),
           x = "Patient ID",
           y = "Expression") +
      ylim(min(df$expression), max(df$expression))
    
    # Print the plot to the PDF
    print(plot)
    
  }, error = function(e) {
    message("Error plotting gene: ", gene, " - ", e$message)
  })
}


grDevices::graphics.off()

# 
# # Iterate over each gene
# for (gene in genes_of_interest) {
#   tryCatch({
#     # Find the correct row for that gene in the count matrix
#     gene_index <- which(rownames(count_mat) == gene)
#     
#     # Extract expression for that gene
#     expression_value <- count_mat[gene_index, ]
#     
#     df <- data.frame(
#       expression = expression_value,
#       Pt_ID = meta.data$Pt_ID,
#       ADpath = meta.data$ADpath
#     )
#     
#     ordering <- c(names(color_mapping)[color_mapping == "red"],
#                   names(color_mapping)[color_mapping == "blue"])
#     df$Pt_ID <- factor(df$Pt_ID, levels = ordering)
#     
#     # Create a data frame of the boxplot statistics
#     boxplot_stats <- df %>%
#       group_by(Pt_ID) %>%
#       summarise(
#         lower_whisker = quantile(expression, 0.25) - 1.5 * IQR(expression),
#         upper_whisker = quantile(expression, 0.75) + 1.5 * IQR(expression)
#       )
#     
#     # Find the maximum value of the upper whiskers
#     max_upper_whisker <- max(boxplot_stats$upper_whisker)
#     
#     k <- which(names(eSVD_obj$teststat_vec) == gene)
#     plot <- ggplot(df, aes(x = Pt_ID, y = expression, fill = Pt_ID)) +
#       ggrastr::rasterize(geom_boxplot(outlier.shape = NA), dpi = 72) + 
#       scale_fill_manual(values = color_mapping) +
#       theme_minimal() +
#       theme(legend.position = "none") +
#       labs(title = paste("Expression of", gene, "\n eSVD LFC: ", round(sereut_obj[k, "logFC_ADpathyes"], 2)), 
#            x = "Donor", 
#            y = "Expression") +
#       ylim(min(df$expression), max_upper_whisker) # Adjust the y-axis limits
#     
#     print(plot)
#   }, error = function(e) {
#     message("Error plotting gene: ", gene, " - ", e$message)
#   })
# }
# 
# # Close PDF device
# grDevices::graphics.off()


#########################
# load("~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_NEBULA.RData")
# pvalue_vec <- nebula_res$summary$p_ADpathyes
# names(pvalue_vec) <- nebula_res$summary$gene
# pvalue_vec_adjusted <- stats::p.adjust(pvalue_vec, method = "BH")
# # Find genes with adjusted p-values <= 0.05
# genes_of_interest <- names(pvalue_vec_adjusted)[pvalue_vec_adjusted <= 0.05]
# genes_of_interest <- intersect(genes_of_interest, rownames(count_mat))
# genes_of_interest

# 
# grDevices::pdf(file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup7/Writeup7_nebula_boxplots.pdf"),
#                width = 8, height = 8, onefile = TRUE)
# # Iterate over each gene
# 
# for (gene in genes_of_interest) {
#   tryCatch({
#     # Find the correct row for that gene in the count matrix
#     gene_index <- which(rownames(count_mat) == gene)
#     
#     # Extract expression for that gene
#     expression_value <- count_mat[gene_index, ]
#     
#     df <- data.frame(
#       expression = expression_value,
#       Pt_ID = meta.data$Pt_ID,
#       ADpath = meta.data$ADpath
#     )
#     
#     ordering <- c(names(color_mapping)[color_mapping == "red"],
#                   names(color_mapping)[color_mapping == "blue"])
#     df$Pt_ID <- factor(df$Pt_ID, levels = ordering)
#     
#     # Create a data frame of the boxplot statistics
#     boxplot_stats <- df %>%
#       group_by(Pt_ID) %>%
#       summarise(
#         lower_whisker = quantile(expression, 0.25) - 1.5 * IQR(expression),
#         upper_whisker = quantile(expression, 0.75) + 1.5 * IQR(expression)
#       )
#     
#     # Find the maximum value of the upper whiskers
#     max_upper_whisker <- max(boxplot_stats$upper_whisker)
#     
#     k <- which(nebula_res$summary$gene == gene)
#     plot <- ggplot(df, aes(x = Pt_ID, y = expression, fill = Pt_ID)) +
#       ggrastr::rasterize(geom_boxplot(outlier.shape = NA), dpi = 72) + 
#       scale_fill_manual(values = color_mapping) +
#       theme_minimal() +
#       theme(legend.position = "none") +
#       labs(title = paste("Expression of", gene, "\nNebula LFC: ", round(nebula_res$summary[k, "logFC_ADpathControl"], 2)), 
#            x = "Donor", 
#            y = "Expression") +
#       ylim(min(df$expression), max_upper_whisker) # Adjust the y-axis limits
#     
#     print(plot)
#   }, error = function(e) {
#     message("Error plotting gene: ", gene, " - ", e$message)
#   })
# }
# 
# # Close PDF device
# grDevices::graphics.off()
