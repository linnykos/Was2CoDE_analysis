# Libraries
rm(list=ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData") 

count_mat <- SeuratObject::LayerData(seurat_obj, layer = "data", assay = "RNA") # genes by cells

# first, create a list (one element per donor) that contains the indicies of all the cells for that donor
meta_data <- seurat_obj@meta.data
donor_ids <- meta_data$donor_id
cell_names <- colnames(count_mat)
cell_df <- data.frame(cell_name = cell_names, donor_id = donor_ids, stringsAsFactors = FALSE)

# Remove "Reference" donors using subset function
seurat_obj <- subset(seurat_obj, subset = ADNC != "Reference")

# Reclassify ADNC levels
seurat_obj$ADNC <- with(seurat_obj@meta.data, 
                        ifelse(ADNC %in% c("NotAD", "Low"), "Control", 
                               ifelse(ADNC %in% c("Intermediate", "High"), "Case", ADNC)))

# then, assign colors to each donor (red = AD, blue = control)

color_mapping <- seurat_obj@meta.data %>%
  select(donor_id, ADNC) %>%
  distinct() %>%
  mutate(color = ifelse(ADNC == "Case", "red", "blue"))

# named vector
color_mapping <- setNames(color_mapping$color, color_mapping$donor_id)

#####
load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_NEBULA.RData")
pvalue_vec <- nebula_res$summary$p_ADNCControl
names(pvalue_vec) <- nebula_res$summary$gene
pvalue_vec_adjusted <- stats::p.adjust(pvalue_vec, method = "BH")
# Find genes with adjusted p-values <= 0.05
genes_of_interest <- names(pvalue_vec_adjusted)[pvalue_vec_adjusted <= 0.05]
genes_of_interest <- intersect(genes_of_interest, rownames(count_mat))

genes_of_interest


grDevices::pdf(file = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup6/Writeup6_nebula_boxplots.pdf"),
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
      donor_id = meta.data$donor_id,
      ADNC = meta.data$ADNC
    )
    
    ordering <- c(names(color_mapping)[color_mapping == "red"],
                  names(color_mapping)[color_mapping == "blue"])
    df$donor_id <- factor(df$donor_id, levels = ordering)
    
    k <- which(nebula_res$summary$gene == gene)
    plot <- ggplot(df, aes(x = donor_id, y = expression, fill = donor_id)) +
      geom_boxplot() + 
      scale_fill_manual(values = color_mapping) +
      theme_minimal() +
      labs(title = paste("Expression of", gene, "\nNebula LFC: ", round(nebula_res$summary[k, "logFC_ADNCControl"], 2)), 
           x = "Donor", 
           y = "Expression")
    
    print(plot)
  }, error = function(e) {
    message("Error plotting gene: ", gene, " - ", e$message)
  })
}

# Close PDF device
grDevices::graphics.off()