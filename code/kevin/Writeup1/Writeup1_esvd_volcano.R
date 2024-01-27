rm(list=ls())
library(EnhancedVolcano)
library(openxlsx)
set.seed(10)

load("../../../../../out/kevin/Writeup1/Writeup1_esvd.RData")

quantile(eSVD_obj$pvalue_list$fdr_vec)
quantile(eSVD_obj$pvalue_list$log10pvalue)

lfc_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
pvalue_vec <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
pval_adj_vec <- eSVD_obj$pvalue_list$fdr_vec
idx <- which(pval_adj_vec <= 0.05)
pCutoff <- max(pvalue_vec[idx])
FCcutoff <- quantile(abs(lfc_vec), probs = 0.9)
xlim <- c(-1,1) * quantile(abs(lfc_vec), probs = 0.99)

idx <- which(abs(pvalue_vec) <= max(xlim))
ylim <- c(0, max(-log10(pvalue_vec[idx])))

res <- data.frame(gene = names(lfc_vec),
                  lfc = lfc_vec,
                  pvalue = pvalue_vec)

plot1 <- EnhancedVolcano::EnhancedVolcano(
  res,
  lab = res$gene,
  x = "lfc",
  y = "pvalue",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  xlim = xlim,
  ylim = ylim
)

ggplot2::ggsave(filename = "../../../figures/kevin/Writeup1/Writeup1_esvd_volcano.png",
                plot1, device = "png", width = 7, height = 7, units = "in")

#########################

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "../../../../../data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
bool_vec <- rep(FALSE, nrow(res))
names(bool_vec) <- res$gene
bool_vec[which(res$gene %in% sun_genes)] <- TRUE
table(bool_vec)

# create a custom volcano plot 

lfc_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
xlim <- c(-1,1) * quantile(abs(lfc_vec), probs = 0.99)

pvalue_vec <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
logpval_vec <- -log10(pvalue_vec)
df <- data.frame(lfc = lfc_vec,
                 log10pval = logpval_vec,
                 name = res$gene,
                 labeling = bool_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "FALSE"), which(df[,"labeling"] == "TRUE")),]

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lfc, 
                                          y = log10pval))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "TRUE"),
                                          ggplot2::aes(label = name, color = labeling),
                                          box.padding = ggplot2::unit(0.5, 'lines'),
                                          point.padding = ggplot2::unit(1.6, 'lines'),
                                          max.overlaps = 50)
plot1 <- plot1 + ggplot2::xlim(xlim)
if(any(pvalue_vec <= 0.05)) {
  plot1 <- plot1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                                       color = "red", linewidth=2)
}
plot1 <- plot1 + ggplot2::ggtitle(paste0("eSVD-DE volcano plot")) +
  ggplot2::xlab("Log fold change") + ggplot2::ylab("eSVD-DE p-value (-Log10)")
plot1 <- plot1 + Seurat::NoLegend()

ggplot2::ggsave(filename = "../../../figures/kevin/Writeup1/Writeup1_esvd_volcano-annotated.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

