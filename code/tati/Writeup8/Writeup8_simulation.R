rm(list=ls())
library(SeuratObject)
library(Seurat)
library(IdeasCustom)
library(devtools)
library(data.table)
library(foreach)
library(doRNG)
library(doParallel)
library(ggplot2)
library(dplyr)
library(reshape2)
registerDoParallel(cores = 4)

set.seed(10)

n_donors <- 20 
n_cells_per_donor <- 200
n_genes <- 3

# # Mean and standard deviation vectors for cases and controls, change mean/sd for other settings
# case_mean_vec <- 0 + runif(n_donors/2, min = -0.5, 0.5)
# case_sd_vec <- rep(5, n_donors/2)
# case_size_vec <- rep(1, n_donors/2)
# control_mean_vec <- 0 + runif(n_donors/2, min = -0.5, 0.5)
# control_sd_vec <- rep(1, n_donors/2)
# control_size_vec <- rep(1, n_donors/2)
# # to do: documentation for the 4 settings
# 
# # Function to generate count data for a single donor
# generate_donor_data <- function(mean_val, sd_val, size_val, n_cells, n_genes) {
#   x <- matrix(rnorm(n_cells * n_genes, mean = mean_val, sd = sd_val), nrow = n_genes, ncol = n_cells)
#   # x <- x - min(x)  # Shift to non-negative values
#   # x <- x / max(x) * size_val * 1000  # Scale
#   # round(x)  # Round to integer counts
# }
# 
# # Generate data for case and control groups
# case_data <- lapply(1:(n_donors/2), function(i) generate_donor_data(case_mean_vec[i], case_sd_vec[i], case_size_vec[i], n_cells_per_donor, n_genes))
# control_data <- lapply(1:(n_donors/2), function(i) generate_donor_data(control_mean_vec[i], control_sd_vec[i], control_size_vec[i], n_cells_per_donor, n_genes))

##########################################################
###different simulation settings####
##########################################################

# 1. Group level: Case vs Control differences
case_group_mean <- 2  # Higher expression in case group
control_group_mean <- 0

# 2. Donor level: Variability between donors within each group
donor_sd <- 0.5

# 3. Cell level: Variability between cells within each donor
cell_sd <- 1

# Generate mean expression levels for each donor
case_mean_vec <- rnorm(n_donors/2, mean = case_group_mean, sd = donor_sd)
control_mean_vec <- rnorm(n_donors/2, mean = control_group_mean, sd = donor_sd)

# Standard deviation vectors for cell-to-cell variability
case_sd_vec <- rep(cell_sd, n_donors/2)
control_sd_vec <- rep(cell_sd, n_donors/2)

# Size factors (if needed for count data generation)
case_size_vec <- rep(1, n_donors/2)
control_size_vec <- rep(1, n_donors/2)

# Documentation for the 4 settings
# 1. Group level difference:
#    - Cases have a mean expression of 2, controls have a mean of 0
# 
# 2. Donor-to-donor variability:
#    - Within each group, donor means are drawn from a normal distribution
#    - The standard deviation (donor_sd = 0.5) determines how much donors vary within their group
# 
# 3. Cell-to-cell variability:
#    - Within each donor, cell expression values are drawn from a normal distribution
#    - (cell_sd = 1) how much cells vary within a donor
# 
# 4. Size factors:
#    - Currently set to 1 for all donors, assuming equal sequencing depth
#    - Can be adjusted if you want to simulate varying sequencing depths between donors

generate_donor_data <- function(donor_mean, cell_sd, size_val, n_cells, n_genes) {
  x <- matrix(rnorm(n_cells * n_genes, mean = donor_mean, sd = cell_sd), 
              nrow = n_genes, 
              ncol = n_cells)
  # Optional: transform to count data if needed
  # x <- x - min(x)  # Shift to non-negative values
  # x <- x / max(x) * size_val * 1000  # Scale
  # x <- round(x)  # Round to integer counts
  return(x)
}

case_data <- lapply(1:(n_donors/2), function(i) 
  generate_donor_data(case_mean_vec[i], case_sd_vec[i], case_size_vec[i], n_cells_per_donor, n_genes))
control_data <- lapply(1:(n_donors/2), function(i) 
  generate_donor_data(control_mean_vec[i], control_sd_vec[i], control_size_vec[i], n_cells_per_donor, n_genes))

########################################################

# Combine all data
count_matrix <- do.call(cbind, c(case_data, control_data))

meta_data <- data.frame(
  donor_id = rep(1:n_donors, each = n_cells_per_donor),
  ADNC = rep(c(rep(1, n_donors/2), rep(0, n_donors/2)), each = n_cells_per_donor),
  nCount_RNA = colSums(count_matrix)
)

# Add row and column names
rownames(count_matrix) <- paste0("Gene_", 1:n_genes)
colnames(count_matrix) <- paste0("Cell_", 1:(n_donors * n_cells_per_donor))
rownames(meta_data) <- colnames(count_matrix)

# Create Seurat object
# seurat_obj <- CreateSeuratObject(counts = count_matrix, meta.data = meta_data)

# Verify the Seurat object
# print(seurat_obj)
# print(table(seurat_obj$ADNC, seurat_obj$donor_id))

# Perform standard Seurat preprocessing
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)

# Extract the processed data matrix
# count_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data"))

# Create meta_cell dataframe
meta_cell <- data.frame(
  cell_id = colnames(count_matrix),
  individual = meta_data$donor_id,
  donor_id = meta_data$donor_id,
  nCount_RNA = meta_data$nCount_RNA,
  ADNC = meta_data$ADNC
)

# Create meta_ind dataframe
meta_ind <- unique(data.frame(
  "individual" = meta_data$donor_id,
  "ADNC" = meta_data$ADNC
))
 
# Handle missing data in meta_ind
for(j in 1:ncol(meta_ind)){
  if(!is.numeric(meta_ind[,j])) next()
  meta_ind[which(is.na(meta_ind[,j])), j] <- median(meta_ind[,j], na.rm = TRUE)
}

# Variables for testing
var2test <- "ADNC"
var2adjust <- setdiff(colnames(meta_ind), c("individual", "ADNC"))
var2test_type <- "binary"
var_per_cell <- "nCount_RNA"

print(head(meta_cell))
print(meta_ind)

# Save intermediate data
save(meta_cell, meta_ind,
     var_per_cell, var2test, var2test_type, var2adjust,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_simulated_data.RData")
   

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- "Was2 analysis on simulated microglia data."

save(meta_cell,
     meta_ind,
     dist_list,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_simulated_microglia_results.RData")

print("Done! :)")

#########################################
library(ggplot2)
library(dplyr)
library(reshape2)
load("~/kzlinlab/projects/subject-de/out/tati/Writeup8/Writeup8_simulated_microglia_results.RData")
ls()

location_data <- dist_list$location
# location_data <- dist_list$distance
# location_data <- dist_list$size

# Prepare data for plotting
plot_data <- data.frame()
for (i in 1:dim(location_data)[1]) {  # For each gene
  gene_data <- location_data[i,,]
  gene_data_melted <- melt(gene_data)
  gene_data_melted$Gene <- dimnames(location_data)[[1]][i]
  
  vec <- sapply(1:nrow(gene_data_melted), function(i){
    x <- gene_data_melted[i,]
    if(x$Var1 >= x$Var2) return("Ignore") ## We might change IdeasCustom::ideas_dist_custom so that we only compute this once (not twice)
    if(x$Var1 <= 10 & x$Var2 <= 10){
      return("Case-Case")
    } else if (x$Var1 > 10 & x$Var2 > 10){
      return("Control-Control")
    } else{
      return("Case-Control")
    }
  })
  
  gene_data_melted$Condition <- vec
  plot_data <- rbind(plot_data, gene_data_melted)
}
plot_data <- plot_data[plot_data$Condition != "Ignore",]

# Rename columns
colnames(plot_data) <- c("Individual1", "Individual2", "Location", "Gene", "Condition")

# Calculate summary statistics
summary_data <- plot_data %>%
  group_by(Gene, Condition) %>%
  summarise(
    Mean = mean(Location),
    SD = sd(Location)
  )

# color palettes
n <- length(unique(plot_data$Gene))
dementia_color_palette <- colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                             rgb(244, 84, 84, maxColorValue = 255)))(n)
no_dementia_color_palette <- colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                rgb(27, 198, 245, maxColorValue = 255)))(n)

p <- ggplot(summary_data, aes(x = Gene, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("Case-Case" = "red", 
                               "Case-Control" = "blue",
                               "Control-Control" = "yellow")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(x = "Gene", y = "Mean Location", title = "Gene Expression Comparison: Dementia vs No Dementia")

print(p)

ggsave("simulation_plot_updated.png", p, width = 12, height = 8, dpi = 300)

# Density plot for a specific gene (e.g., Gene-1)
gene_of_interest <- "Gene_1"
gene_data <- plot_data %>% filter(Gene == gene_of_interest)

p_density <- ggplot(gene_data, aes(x = Location, fill = Condition)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("Case-Case" = "red", 
                               "Case-Control" = "blue",
                               "Control-Control" = "yellow")) +
  theme_minimal() +
  # ylim(c(0,5)) + # (Include this if you need to zoom in)
  labs(x = "Location", y = "Density", title = paste("Distribution of", gene_of_interest, "Expression"))

# Display the density plot
print(p_density)

ggsave(paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8/Writeup8_", "simulation_boxplot_updated.png"), p, width = 12, height = 8, dpi = 300)
ggsave(paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup8/Writeup8_", "density_plot_gene1.png"), p_density, width = 10, height = 6, dpi = 300)


