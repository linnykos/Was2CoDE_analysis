# Libraries
rm(list=ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

count_mat <- SeuratObject::LayerData(ss_data_norm, layer = "data", assay = "RNA") # genes by cells

# first, create a list (one element per donor) that contains the indicies of all the cells for that donor
meta_data <- ss_data_norm@meta.data
Pt_ids <- meta_data$Pt_ID
cell_names <- colnames(count_mat)
cell_df <- data.frame(cell_name = cell_names, Pt_ID = Pt_ids, stringsAsFactors = FALSE)


# then, assign colors to each donor (red = AD, blue = control)

color_mapping <- meta_data %>%
  select(Pt_ID, Study_Designation) %>%
  distinct() %>%
  mutate(color = ifelse(Study_Designation == "AD", "red", "blue"))

# named vector
color_mapping <- setNames(color_mapping$color, color_mapping$Pt_ID)

#####

load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_NEBULA.RData")
pvalue_vec <- nebula_res$summary$p_Study_DesignationAD
names(pvalue_vec) <- nebula_res$summary$gene
pvalue_vec_adjusted <- stats::p.adjust(pvalue_vec, method = "BH")
# Find genes with adjusted p-values <= 0.05
genes_of_interest <- names(pvalue_vec_adjusted)[pvalue_vec_adjusted <= 0.05]
genes_of_interest <- intersect(genes_of_interest, rownames(count_mat))

genes_of_interest


grDevices::pdf(file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/Writeup5_nebula_boxplots.pdf"),
    width = 8, height = 8, onefile = TRUE)
# Iterate over each gene

  for (gene in genes_of_interest) {
    tryCatch({
      # Find the correct row for that gene in the count matrix
      gene_index <- which(rownames(count_mat) == gene)
      
      # Extract expression for that gene
      expression_value <- count_mat[gene_index, ]
      
      df <- data.frame(
        expression = expression_value,
        Pt_ID = meta_data$Pt_ID,
        Study_Designation = meta_data$Study_Designation
      )
      
      ordering <- c(names(color_mapping)[color_mapping == "red"],
                    names(color_mapping)[color_mapping == "blue"])
      df$Pt_ID <- factor(df$Pt_ID, levels = ordering)
      
      k <- which(nebula_res$summary$gene == gene)
      plot <- ggplot(df, aes(x = Pt_ID, y = expression, fill = Pt_ID)) +
        geom_boxplot() + 
        scale_fill_manual(values = color_mapping) +
        theme_minimal() +
        labs(title = paste("Expression of", gene, "\nNebula LFC: ", round(nebula_res$summary[k, "logFC_Study_DesignationAD"], 2)), 
             x = "Donor", 
             y = "Expression")
      
      print(plot)
    }, error = function(e) {
      message("Error plotting gene: ", gene, " - ", e$message)
    })
  }
  
  # Close PDF device
  grDevices::graphics.off()