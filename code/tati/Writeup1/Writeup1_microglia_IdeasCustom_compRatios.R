library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideascustom.RData")

gene_index <- 738  #gene index for "SMO"

# Extract gene data from the dist_list
gene_distance <- dist_list$distance[gene_index, ,]
gene_location <- dist_list$location[gene_index, ,]
gene_size <- dist_list$size[gene_index, ,]
gene_shape <- dist_list$shape[gene_index, ,] 

# Initialize matrices to store results
n <- nrow(meta_ind)  # Number of donors
location_ratio <- matrix(NA, n, n)
size_ratio <- matrix(NA, n, n)
shape_ratio <- matrix(NA, n, n)

# Calculate the ratios for each pair of donors
for (i in 1:n) {
  for (j in 1:n) {
    distance_squared <- gene_distance[i, j]^2
      location_diff <- gene_location[i,j]
      size_diff <- gene_size[i,j]
      shape_diff <- gene_shape[i,j] 
      
      location_ratio[i, j] <- location_diff / distance_squared
      size_ratio[i, j] <- size_diff / distance_squared
      shape_ratio[i, j] <- shape_ratio[i, j] <- (distance_squared - location_diff - size_diff) / distance_squared
  }
}
rownames(location_ratio) <- meta_ind$individual
rownames(size_ratio) <- meta_ind$individual
rownames(shape_ratio) <- meta_ind$individual
colnames(location_ratio) <- meta_ind$individual
colnames(size_ratio) <- meta_ind$individual
colnames(shape_ratio) <- meta_ind$individual

dementia_indices <- which(meta_ind$CognitiveStatus == "Dementia")
no_dementia_indices <- which(meta_ind$CognitiveStatus == "No_dementia")

location_ratio_dnd <- location_ratio[dementia_indices, no_dementia_indices]
size_ratio_dnd <- size_ratio[dementia_indices, no_dementia_indices]
shape_ratio_dnd <- shape_ratio[dementia_indices, no_dementia_indices]

location_ratio_dd <- location_ratio[dementia_indices, dementia_indices]
size_ratio_dd <- size_ratio[dementia_indices, dementia_indices]
shape_ratio_dd <- shape_ratio[dementia_indices, dementia_indices]

location_ratio_ndnd <- location_ratio[no_dementia_indices, no_dementia_indices]
size_ratio_ndnd <- size_ratio[no_dementia_indices, no_dementia_indices]
shape_ratio_ndnd <- shape_ratio[no_dementia_indices, no_dementia_indices]

location_data_dnd <- as.data.frame(location_ratio_dnd) %>% mutate(Category = "Case & Control",RatioType = "Location")
location_data_dd <- as.data.frame(location_ratio_dd) %>% mutate(Category = "Case & Case", RatioType = "Location")
location_data_ndnd <- as.data.frame(location_ratio_ndnd) %>% mutate(Category = "Control & Control", RatioType = "Location")
size_data_dnd <- as.data.frame(size_ratio_dnd) %>% mutate(Category = "Case & Control", RatioType = "Size")
size_data_dd <- as.data.frame(size_ratio_dd) %>% mutate(Category = "Case & Case", RatioType = "Size")
size_data_ndnd <- as.data.frame(size_ratio_ndnd) %>% mutate(Category = "Control & Control", RatioType = "Size")
shape_data_dnd <- as.data.frame(shape_ratio_dnd) %>% mutate(Category = "Case & Control", RatioType = "Shape")
shape_data_dd <- as.data.frame(shape_ratio_dd) %>% mutate(Category = "Case & Case", RatioType = "Shape")
shape_data_ndnd <- as.data.frame(shape_ratio_ndnd) %>% mutate(Category = "Control & Control", RatioType = "Shape")

location_combined_data <- bind_rows(location_data_dnd, location_data_dd, location_data_ndnd)
size_combined_data <- bind_rows(size_data_dnd,size_data_dd, size_data_ndnd)
shape_combined_data <- bind_rows(shape_data_dnd, shape_data_dd,shape_data_ndnd)

location_long_data <- location_combined_data %>%
  pivot_longer(cols = -c(Category, RatioType), names_to = "Pair", values_to = "LocationRatio")
size_long_data <- size_combined_data %>%
  pivot_longer(cols = -c(Category, RatioType), names_to = "Pair", values_to = "SizeRatio")
shape_long_data <- shape_combined_data %>%
  pivot_longer(cols = -c(Category, RatioType), names_to = "Pair", values_to = "ShapeRatio")


plot_location <- ggplot(location_long_data, aes(x = LocationRatio, fill = Category)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Case & Control" = "blue", "Case & Case" = "red", "Control & Control" = "green")) +
  labs(title = "Distribution of Location/Distance^2 Ratios", x = "Location/Distance^2 Ratio", y = "Density") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_size <- ggplot(size_long_data, aes(x = SizeRatio, fill = Category)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Case & Control" = "blue", "Case & Case" = "red", "Control & Control" = "green")) +
  labs(title = "Distribution of Size Ratios", x = "Size Ratio", y = "Density") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_shape <- ggplot(shape_long_data, aes(x = ShapeRatio, fill = Category)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Case & Control" = "blue", "Case & Case" = "red", "Control & Control" = "green")) +
  labs(title = "Distribution of Shape Ratios", x = "Shape Ratio", y = "Density") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-histogram_plot_location_IdeasCustom_1.png"),
                plot_location, device = "png", width = 7, height = 7, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-histogram_plot_size_IdeasCustom_1.png"),
                plot_size, device = "png", width = 7, height = 7, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-histogram_plot_shape_IdeasCustom_1.png"),
                plot_shape, device = "png", width = 7, height = 7, units = "in")


boxplot_location <-
  ggplot(location_long_data, aes(x=Category, y=LocationRatio, fill=Category)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Boxplot of Location Ratios") +
  xlab("Category") +
  ylab("Location Ratio")


boxplot_size<-
  ggplot(size_long_data, aes(x=Category, y=SizeRatio, fill=Category)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Boxplot of Size Ratios") +
  xlab("Category") +
  ylab("Size Ratio")

boxplot_shape<-
ggplot(shape_long_data, aes(x=Category, y=ShapeRatio, fill=Category)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Boxplot of Shape Ratios") +
  xlab("Category") +
  ylab("Shape Ratio")

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-histogram_boxplot_location_IdeasCustom_1.png"),
                boxplot_location, device = "png", width = 5, height = 7, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-histogram_boxplot_size_IdeasCustom_1.png"),
                boxplot_size, device = "png", width = 5, height = 7, units = "in")
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de-histogram_boxplot_shape_IdeasCustom_1.png"),
                boxplot_shape, device = "png", width = 5, height = 7, units = "in")