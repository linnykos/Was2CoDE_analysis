# Libraries
rm(list=ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)


load("~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData")
count_mat <- SeuratObject::LayerData(ss_data_norm, layer = "data", assay = "RNA") # genes by cells

# first, create a list (one element per donor) that contains the indicies of all the cells for that donor
meta_data <- ss_data_norm@meta.data
Pt_ids <- meta_data$Pt_ID
cell_names <- colnames(count_mat)
cell_df <- data.frame(cell_name = cell_names, Pt_ID = Pt_ids, stringsAsFactors = FALSE)


# then, assign colors to each donor (red = AD, blue = control)

color_mapping <- meta_data %>%
  select(Pt_ID, CognitiveStatus) %>%
  distinct() %>%
  mutate(color = ifelse(CognitiveStatus == "Dementia", "red", "blue"))

# named vector
color_mapping <- setNames(color_mapping$color, color_mapping$Pt_ID)

#####

load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_nebula.RData")
pvalue_vec <- nebula_res$summary[,"p_CognitiveStatusNo dementia"]
names(pvalue_vec) <- nebula_res$summary[,"gene"]
pvalue_vec <- stats::p.adjust(pvalue_vec, method = "BH")
genes_of_interest <- names(pvalue_vec)[pvalue_vec <= 0.05]
genes_of_interest <- intersect(genes_of_interest, rownames(count_mat))

grDevices::pdf(file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_nebula_boxplots.pdf"),
    width = 8, height = 8, onefile = TRUE)
# Iterate over each gene
for (gene in genes_of_interest) {
  # Find the correct row for that gene in the count matrix
  gene_index <- which(rownames(count_mat) == gene)
  
  # Extract expression for that gene
  expression_value <- count_mat[gene_index, ]

  df <- data.frame(
    expression = expression_value,
    Pt_ID = meta_data$Pt_ID,
    CognitiveStatus = meta_data$CognitiveStatus
  )
  
  ordering <- c(names(color_mapping)[color_mapping == "red"],
                names(color_mapping)[color_mapping == "blue"])
  df$Pt_ID <- factor(df$Pt_ID, levels = ordering)
  
  mapping <- c(Dementia = "red", No_dementia = "blue")
  k <- which(nebula_res$summary$gene == gene)
  plot <- ggplot(df, aes(x = Pt_ID, y = expression, fill = Pt_ID)) +
    geom_boxplot() + 
    # ggrastr::rasterise(geom_boxplot(), dpi = 120) +
    scale_fill_manual(values = color_mapping) +
    theme_minimal() +
    labs(title = paste("Expression of", gene, "\nNebula LFC: ", round(nebula_res$summary[k, "logFC_CognitiveStatusNo dementia"], 2)), x = "Donor", y = "Expression")
  
  print(plot)
}

grDevices::graphics.off()

