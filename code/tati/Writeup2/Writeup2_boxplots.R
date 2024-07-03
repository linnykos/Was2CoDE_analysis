# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_nebula.RData")

# Extract the summary data frame
nebula_summary <- nebula_res$summary

# Define the significant genes from the volcano plot (example genes, replace with actual genes from your plot)
significant_genes <- c("SORCS1", "DPYS", "ASTN1", "GLT1D1", "BRINP2", "PBX1", "MAPK10", "RP11-762E8.1", "FSTL5", "ARHGAP10")

# Filter the summary data frame to include only significant genes
filtered_data <- nebula_summary %>%
  filter(gene %in% significant_genes)

# Create a boxplot for the significant genes
plot1 <- ggplot(filtered_data, aes(x = gene, y = `logFC_(Intercept)`, fill = gene)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ggtitle("Boxplot of Significant Genes Nebula") +
  xlab("Genes") +
  ylab("logFC (Intercept)")


ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_nebula_boxplot.png"),
                plot1, device = "png", width = 5, height = 7, units = "in")
