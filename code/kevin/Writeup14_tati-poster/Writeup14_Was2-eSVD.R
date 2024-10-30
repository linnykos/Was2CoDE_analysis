rm(list=ls())

library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)

plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup14/"
df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Prater_dataset_ingredients.csv")
rownames(df) <- df$Gene

x_pval <- df[,"eSVD_pval"]
names(x_pval) <- rownames(df)
y_pval <- df[,"Was2_pval"]
names(y_pval) <- rownames(df)
cor(-log10(x_pval), -log10(y_pval))

x_padj <- stats::p.adjust(x_pval, method = "BH")
y_padj <- stats::p.adjust(y_pval, method = "BH")
x_genes <- names(x_padj)[x_padj <= 0.05]
y_genes <- names(y_padj)[y_padj <= 0.1]

x_fdrcutoff <- -log10(max(x_pval[x_genes]))
y_fdrcutoff <- -log10(max(y_pval[y_genes]))

ggplot_df <- data.frame(
  Gene = rownames(df),
  x = -log10(x_pval),
  y = -log10(y_pval)
)

label_bool_vec <- rep("", nrow(ggplot_df))
names(label_bool_vec) <- rownames(ggplot_df) 
any_genes <- unique(c(x_genes, y_genes))
label_bool_vec[any_genes] <- any_genes
ggplot_df$label <- label_bool_vec
ggplot_df$color <- sapply(label_bool_vec, function(x){
  nchar(x) > 0
})

color_vec <- setNames(c("gray30", "#bc6a17"), 
                      c(FALSE, TRUE))

# make the plot
plot1 <- ggplot2::ggplot(data = ggplot_df, mapping = aes(x = x,
                                                         y = y,
                                                         label = label,
                                                         color = color)) + 
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray") +    
  geom_vline(xintercept = 0, linewidth = 0.5, color = "gray") + 
  scale_color_manual(values = color_vec) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = y_fdrcutoff, linetype = "dashed") +  
  geom_vline(xintercept = x_fdrcutoff, linetype = "dashed") +
  geom_text_repel(size = 2, colour = "black") +
  labs(x = paste0("eSVD-DE (Signed -log10 p-value)"), 
       y = paste0("Was2 (Signed -log10 p-value)")) + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),        # Smaller title text size
    axis.title.x = element_text(size = 10),      # Smaller x-axis label text size
    axis.title.y = element_text(size = 10)       # Smaller y-axis label text size
  ) +
  Seurat::NoLegend()

ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_Was2-eSVD.png"),
       height = 3.1,
       width = 3)