rm(list=ls())

df <- read.csv("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Prater_dataset_ingredients.csv")
rownames(df) <- df$Gene
load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_microglia_ideascustom.RData")

names(dist_list)
dim(dist_list[[1]])
dist_list[[1]][1,1:3,1:3]
dimnames(dist_list[[1]])

# grab the %mean for all the genes in Wasserstein

# first compute the average mean/was2 for each person
tri_mask <- upper.tri(dist_list[[1]][1,,], diag = FALSE)
gene_perc_mat <- sapply(1:dim(dist_list[[1]])[1], function(j){
  was2_mat <- dist_list[["distance"]][j,,]^2
  loc_mat <- dist_list[["location"]][j,,]
  perc_mat <- loc_mat/was2_mat
  perc_mat[tri_mask]
})
gene_perc_vec <- Matrix::colMeans(gene_perc_mat)
names(gene_perc_vec) <- dimnames(dist_list[[1]])[[1]]

gene_intersection <- intersect(names(gene_perc_vec),
                               df$Gene)
deseq2_logfc <- abs(df[gene_intersection, "DESeq2_logFC"])
gene_perc_vec <- gene_perc_vec[gene_intersection]

ggplot_df <- data.frame(
  Gene = gene_intersection,
  deseq2 = deseq2_logfc,
  was_perc = gene_perc_vec
)

# make the plot
plot_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup14/"

plot1 <- ggplot2::ggplot(data = ggplot_df, mapping = aes(x = deseq2,
                                                         y = was_perc)) + 
  geom_point() +
  ggtitle(paste0("Cor: ", round(stats::cor(deseq2_logfc, gene_perc_vec), 2))) + 
  theme_minimal() 

ggsave(plot1, 
       filename = paste0(plot_folder, "Writeup14_was2_tmp.png"),
       height = 5,
       width = 5)
  