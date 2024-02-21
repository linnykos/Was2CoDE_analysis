rm(list=ls())
library(EnhancedVolcano)
library(openxlsx)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")

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

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_esvd_volcano.png",
                plot1, device = "png", width = 7, height = 7, units = "in")

#########################

lfc_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
xlim <- c(-1,1) * quantile(abs(lfc_vec), probs = 0.99)

pvalue_vec <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
logpval_vec <- -log10(pvalue_vec)
fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
fdr_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
fdr_bool <- rep(FALSE, length(fdr_vec))
names(fdr_bool) <- names(fdr_vec)
fdr_bool[which(fdr_vec <= 0.05)] <- TRUE

#########################

prater_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/43587_2023_424_MOESM3_ESM.xlsx",
  sheet = "AD_vs_Ctrl_DEGs",
  startRow = 3
) 
prater_genes <- sort(unique(prater_sheet[which(prater_sheet[,"padj"] <= 0.05),"Gene"]))
bool_vec <- rep(FALSE, length(fdr_bool))
names(bool_vec) <- names(fdr_bool)
bool_vec[which(names(bool_vec) %in% prater_genes)] <- TRUE
bool_vec1 <- bool_vec*fdr_bool # both in validation dataset and our FDR
bool_vec2 <- bool_vec - bool_vec1 # only in validation
label_vec <- rep("None", length(bool_vec))
label_vec[which(bool_vec2 == 1)] <- "Validation"
label_vec[which(bool_vec1 == 1)] <- "Both"

# create a custom volcano plot 
df <- data.frame(lfc = lfc_vec,
                 log10pval = logpval_vec,
                 name = res$gene,
                 labeling = factor(label_vec))
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "None"), 
           which(df[,"labeling"] == "Validation"),
           which(df[,"labeling"] == "Both")),]

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lfc, 
                                          y = log10pval))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c(Both = "red",
                                                       None = "black", 
                                                       Validation = "cadetblue2"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "Both"),
                                          ggplot2::aes(label = name, color = labeling),
                                          box.padding = ggplot2::unit(0.5, 'lines'),
                                          point.padding = ggplot2::unit(1.6, 'lines'),
                                          max.overlaps = 50)
plot1 <- plot1 + ggplot2::xlim(xlim)
if(any(pvalue_vec <= 0.05)) {
  plot1 <- plot1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(fdr_vec <= 0.05)]), linetype="dashed", 
                                       color = "red", linewidth=2)
}
plot1 <- plot1 + ggplot2::ggtitle(paste0("eSVD-DE volcano plot, Prater genes")) +
  ggplot2::xlab("Log fold change") + ggplot2::ylab("eSVD-DE p-value (-Log10)")
plot1 <- plot1 + Seurat::NoLegend()

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_esvd_volcano_annotated-prater.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

#########################

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
bool_vec <- rep(FALSE, length(fdr_bool))
names(bool_vec) <- names(fdr_bool)
bool_vec[which(names(bool_vec) %in% sun_genes)] <- TRUE
bool_vec1 <- bool_vec*fdr_bool # both in validation dataset and our FDR
bool_vec2 <- bool_vec - bool_vec1 # only in validation
label_vec <- rep("None", length(bool_vec))
label_vec[which(bool_vec2 == 1)] <- "Validation"
label_vec[which(bool_vec1 == 1)] <- "Both"

# create a custom volcano plot 
df <- data.frame(lfc = lfc_vec,
                 log10pval = logpval_vec,
                 name = res$gene,
                 labeling = factor(label_vec))
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "None"), 
           which(df[,"labeling"] == "Validation"),
           which(df[,"labeling"] == "Both")),]

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lfc, 
                                          y = log10pval))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c(Both = "red",
                                                       None = "black", 
                                                       Validation = "cadetblue2"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "Both"),
                                          ggplot2::aes(label = name, color = labeling),
                                          box.padding = ggplot2::unit(0.5, 'lines'),
                                          point.padding = ggplot2::unit(1.6, 'lines'),
                                          max.overlaps = 50)
plot1 <- plot1 + ggplot2::xlim(xlim)
if(any(pvalue_vec <= 0.05)) {
  plot1 <- plot1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(fdr_vec <= 0.05)]), linetype="dashed", 
                                       color = "red", linewidth=2)
}
plot1 <- plot1 + ggplot2::ggtitle(paste0("eSVD-DE volcano plot, Sun genes")) +
  ggplot2::xlab("Log fold change") + ggplot2::ylab("eSVD-DE p-value (-Log10)")
plot1 <- plot1 + Seurat::NoLegend()

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_esvd_volcano_annotated-sun.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###############

df[intersect(intersect(sun_genes, prater_genes), fdr_genes),]
df["APOE",]
          
