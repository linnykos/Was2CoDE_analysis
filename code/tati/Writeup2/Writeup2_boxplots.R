# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_nebula.RData")
nebula_summary <- nebula_res$summary

significant_genes <- c("SORCS1", "DPYS", "ASTN1", "GLT1D1", "BRINP2", "PBX1", "MAPK10", "RP11-762E8.1", "FSTL5", "ARHGAP10")

filtered_data <- nebula_summary %>%
  filter(gene %in% significant_genes)

for (gene in significant_genes) {
  # Filter the data for the current gene
  gene_data <- nebula_summary %>%
    filter(gene == !!gene)
  
  plot1 <- ggplot(gene_data, aes(x = CognitiveStatus, y = `logFC_CognitiveStatusNo dementia`, fill = CognitiveStatus)) +
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
    xlab("Cognitive Status") +
    ylab("logFC (AD vs. non-AD)")

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup2/Writeup2_nebula_boxplot.png"),
                plot1, device = "png", width = 5, height = 7, units = "in")
}