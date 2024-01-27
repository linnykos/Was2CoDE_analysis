rm(list=ls())
library(EnhancedVolcano)
library(openxlsx)
set.seed(10)

load("../../../../../out/kevin/Writeup1/Writeup1_esvd.RData")
pvalue_esvd_vec <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
logpvalue_esvd_vec <- eSVD_obj$pvalue_list$log10pvalue
pvalue_adj_esvd_vec <- eSVD_obj$pvalue_list$fdr_vec

load("../../../../../out/kevin/Writeup1/Writeup1_nebula.RData")
res <- nebula_res$summary
pvalue_nebula_vec <- res[,"p_CognitiveStatusNo dementia"]
names(pvalue_nebula_vec) <- res$gene
logpvalue_nebula_vec <- -log10(pvalue_nebula_vec)
pvalue_adj_nebula_vec <- stats::p.adjust(pvalue_nebula_vec, method = "BH")

gene_vec <- sort(intersect(names(logpvalue_esvd_vec), names(logpvalue_nebula_vec)))
pvalue_esvd_vec <- pvalue_esvd_vec[gene_vec]
pvalue_adj_esvd_vec <- pvalue_adj_esvd_vec[gene_vec]
logpvalue_esvd_vec <- logpvalue_esvd_vec[gene_vec]
pvalue_nebula_vec <- pvalue_nebula_vec[gene_vec]
pvalue_adj_nebula_vec <- pvalue_adj_nebula_vec[gene_vec]
logpvalue_nebula_vec <- logpvalue_nebula_vec[gene_vec]

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "../../../../../data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
bool_vec <- rep(FALSE, nrow(res))
names(bool_vec) <- gene_vec
bool_vec[which(gene_vec %in% sun_genes)] <- TRUE

df <- data.frame(name = gene_vec,
                 esvd_logp = logpvalue_esvd_vec,
                 nebula_logp = logpvalue_nebula_vec,
                 labeling = bool_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "FALSE"), which(df[,"labeling"] == "TRUE")),]

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = esvd_logp, 
                                          y = nebula_logp))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "TRUE"),
                                          ggplot2::aes(label = name, color = labeling),
                                          box.padding = ggplot2::unit(0.5, 'lines'),
                                          point.padding = ggplot2::unit(1.6, 'lines'),
                                          max.overlaps = 50)
plot1 <- plot1 + ggplot2::xlim(xlim)
if(any(pvalue_adj_esvd_vec <= 0.05)) {
  plot1 <- plot1 + ggplot2::geom_vline(xintercept=min(logpvalue_esvd_vec[which(pvalue_adj_esvd_vec <= 0.05)]), linetype="dashed", 
                                       color = "red", linewidth=2)
}
if(any(pvalue_adj_nebula_vec <= 0.05)) {
  plot1 <- plot1 + ggplot2::geom_hline(yintercept=min(logpvalue_nebula_vec[which(pvalue_adj_nebula_vec <= 0.05)]), linetype="dashed", 
                                       color = "red", linewidth=2)
}
plot1 <- plot1 + ggplot2::ggtitle(paste0("NEBULA vs. eSVD volcano plot\nCorrelation: ", round(stats::cor(logpvalue_esvd_vec, logpvalue_nebula_vec), 2))) +
  ggplot2::xlab("eSVD p-value (-Log10)") + ggplot2::ylab("NEBULA p-value (-Log10)")
plot1 <- plot1 + Seurat::NoLegend()

ggplot2::ggsave(filename = "../../../figures/kevin/Writeup1/Writeup1_nebula-esvd_volcano-annotated.png",
                plot1, device = "png", width = 7, height = 7, units = "in")



