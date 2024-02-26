rm(list=ls())

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_ideas_pvalues.RData")
gene_names <- names(pval_res)
ideas_pvalue <- pval_res
ideas_log10 <- -log10(ideas_pvalue)
ideas_fdr <- stats::p.adjust(ideas_pvalue, method = "BH")
ideas_cutoff <- min(ideas_log10[which(ideas_fdr <= 0.05)])
ideas_selected_genes <- names(ideas_fdr)[which(ideas_fdr <= 0.05)]

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
esvd_pvalue <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
esvd_pvalue <- esvd_pvalue[gene_names]
esvd_log10 <- -log10(esvd_pvalue)
esvd_fdr <- eSVD_obj$pvalue_list$fdr_vec
esvd_cutoff <- min(esvd_log10[which(esvd_fdr <= 0.05)])
esvd_selected_genes <- names(esvd_fdr)[which(esvd_fdr <= 0.05)]

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
sun_genes <- intersect(gene_names, sun_genes)

#################

m <- length(sun_genes)
n <- length(gene_names) - m
k <- length(ideas_selected_genes)
x <- length(intersect(ideas_selected_genes, sun_genes))
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
-log10(fisher)

####

selected_genes <- sort(unique(c(ideas_selected_genes, esvd_selected_genes)))

yellow_col <- rgb(255, 205, 114, maxColorValue = 255) # deseq2
orange_col <- rgb(235, 134, 47, maxColorValue = 255) # esvd
blue_col <- rgb(50, 174, 255, maxColorValue = 255) # ideas
green_col <- rgb(70, 177, 70, maxColorValue = 255) # sun genes

purple_col <- rgb(122, 49, 126, maxColorValue = 255)
red_col <- rgb(255, 90, 90, maxColorValue = 255) 

p <- length(gene_names)
col_vec <- rep(rgb(0.6, 0.6, 0.6), p)
names(col_vec) <- gene_names
col_vec[setdiff(sun_genes, selected_genes)] <- red_col
# col_vec[esvd_selected_genes] <- orange_col
# col_vec[ideas_selected_genes] <- blue_col

size_vec <- rep(1, p)
names(size_vec) <- gene_names
size_vec[setdiff(sun_genes, selected_genes)] <- 0.5

transparent_vec <- rep(TRUE, p)
names(transparent_vec) <- gene_names
transparent_vec[c(selected_genes, sun_genes)] <- FALSE

circle_vec <- rep(FALSE, p)
names(circle_vec) <- gene_names
circle_vec[intersect(sun_genes, selected_genes)] <- TRUE

sun_vec <- rep(FALSE, p)
names(sun_vec) <- gene_names
sun_vec[sun_genes] <- TRUE

df <- data.frame(
  name = gene_names,
  esvd_log10pvalue = esvd_log10,
  ideas_log10pvalue = ideas_log10,
  sun = sun_vec,
  col = col_vec,
  circled = circle_vec,
  transparent = transparent_vec,
  size = size_vec
)
# order_idx <- c(which(df$transparent == T),
#                intersect(which(df$transparent == F), which(!df$name %in% selected_genes)),
#                intersect(which(df$transparent == F), which(df$name %in% selected_genes)))
# df <- df[order_idx,]

# Create the scatterplot using ggplot2
plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = ideas_log10pvalue, 
                                          y = esvd_log10pvalue, 
                                          size = size))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = col, alpha = ifelse(transparent, 0, 1)), show.legend = FALSE)
plot1 <- plot1 + ggplot2::scale_size_continuous(range = c(0.5, 1.5))
plot1 <- plot1 + ggplot2::scale_alpha_continuous(range = c(0.3, 1))
plot1 <- plot1 + ggplot2::scale_color_identity()
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & sun == TRUE), color = red_col, size = 3, shape = 1)
plot1 <- plot1 + ggplot2::geom_hline(yintercept = esvd_cutoff, color = "white", size = 2)
plot1 <- plot1 + ggplot2::geom_vline(xintercept = ideas_cutoff, color = "white", size = 2)
plot1 <- plot1 + ggplot2::geom_hline(yintercept = esvd_cutoff, linetype = "dashed", color = orange_col, size = 1.5)
plot1 <- plot1 + ggplot2::geom_vline(xintercept = ideas_cutoff, linetype = "dashed", color = blue_col, size = 1.5)
plot1 <- plot1 + ggplot2::labs(x = "", y = "", title = "")
plot1 <- plot1 + ggplot2::labs(title = paste0("Correlation: ", round(stats::cor(ideas_log10, esvd_log10), 2)))

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_esvd-ideas_volcano-pretty.png",
                plot1, device = "png", width = 3, height = 4, units = "in")
