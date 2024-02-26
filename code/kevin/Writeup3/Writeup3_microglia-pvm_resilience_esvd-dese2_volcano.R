rm(list=ls())
library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_deseq2.RData")
gene_vec <- rownames(deseq2_res)
deseq2_pvalue <- deseq2_res$pvalue
deseq2_log10 <- -log10(deseq2_pvalue)
deseq2_log10[is.na(deseq2_log10)] <- 0
deseq2_fdr <- deseq2_res$padj
deseq2_cutoff <- min(deseq2_log10[which(deseq2_fdr <= 0.05)])

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
esvd_pvalue <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
esvd_pvalue <- esvd_pvalue[gene_vec]
esvd_log10 <- -log10(esvd_pvalue)
esvd_fdr <- eSVD_obj$pvalue_list$fdr_vec
esvd_cutoff <- min(esvd_log10[which(esvd_fdr <= 0.05)])

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
sun_genes <- intersect(gene_vec, sun_genes)

label_vec <- rep("None", length(gene_vec))
label_vec[which(gene_vec %in% sun_genes)] <- "Validation"
label_vec[intersect(which(deseq2_log10 > deseq2_cutoff),
                    which(gene_vec %in% sun_genes))] <- "Both"
label_vec[intersect(which(esvd_log10 > esvd_cutoff),
                    which(gene_vec %in% sun_genes))] <- "Both"

# cor.test(deseq2_log10, esvd_log10)

# create a custom volcano plot 
df <- data.frame(deseq2_log10 = deseq2_log10,
                 esvd_log10 = esvd_log10,
                 name = gene_vec,
                 labeling = factor(label_vec))
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "None"), 
           which(df[,"labeling"] == "Validation"),
           which(df[,"labeling"] == "Both")),]

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = deseq2_log10, 
                                          y = esvd_log10))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c(Both = "red",
                                                       None = "black", 
                                                       Validation = "cadetblue2"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "Both"),
                                          ggplot2::aes(label = name, color = labeling, wt = 0.5))
plot1 <- plot1 + ggplot2::geom_hline(yintercept=esvd_cutoff, linetype="dashed", 
                                     color = "red", linewidth=2)
plot1 <- plot1 + ggplot2::geom_vline(xintercept=deseq2_cutoff, linetype="dashed", 
                                     color = "red", linewidth=2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Correlation: ", round(stats::cor(deseq2_log10, esvd_log10), 2))) +
  ggplot2::xlab("DESeq2 p-value (-Log10)") + ggplot2::ylab("eSVD-DE p-value (-Log10)")
plot1 <- plot1 + Seurat::NoLegend()

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_esvd-deseq2_volcano.png",
                plot1, device = "png", width = 7, height = 10, units = "in")

####################
####################

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_deseq2.RData")
gene_vec <- rownames(deseq2_res)
deseq2_lfc <- deseq2_res$log2FoldChange
deseq2_lfc[is.na(deseq2_lfc)] <- 0
deseq2_pvalue <- deseq2_res$pvalue
deseq2_log10 <- -log10(deseq2_pvalue)
deseq2_log10[is.na(deseq2_log10)] <- 0
deseq2_fdr <- deseq2_res$padj
deseq2_cutoff <- min(deseq2_log10[which(deseq2_fdr <= 0.05)])

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
esvd_lfc <- eSVD_obj$control_mean - eSVD_obj$case_mean
esvd_lfc <- pmin(esvd_lfc, stats::quantile(esvd_lfc, probs = 0.99))
esvd_lfc <- pmax(esvd_lfc, stats::quantile(esvd_lfc, probs = 0.01))
esvd_pvalue <- 10^(-eSVD_obj$pvalue_list$log10pvalue)
esvd_pvalue <- esvd_pvalue[gene_vec]
esvd_log10 <- -log10(esvd_pvalue)
esvd_fdr <- eSVD_obj$pvalue_list$fdr_vec
esvd_cutoff <- min(esvd_log10[which(esvd_fdr <= 0.05)])

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
sun_genes <- intersect(gene_vec, sun_genes)

label_vec <- rep("None", length(gene_vec))
label_vec[which(gene_vec %in% sun_genes)] <- "Validation"
label_vec[intersect(which(deseq2_log10 > deseq2_cutoff),
                    which(gene_vec %in% sun_genes))] <- "Both"
label_vec[intersect(which(esvd_log10 > esvd_cutoff),
                    which(gene_vec %in% sun_genes))] <- "Both"

# cor.test(deseq2_lfc, esvd_lfc)

# create a custom volcano plot 
df <- data.frame(deseq2_lfc = deseq2_lfc,
                 esvd_lfc = esvd_lfc,
                 name = gene_vec,
                 labeling = factor(label_vec))
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "None"), 
           which(df[,"labeling"] == "Validation"),
           which(df[,"labeling"] == "Both")),]

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = deseq2_lfc, 
                                          y = esvd_lfc))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c(Both = "red",
                                                       None = "black", 
                                                       Validation = "cadetblue2"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "Both"),
                                          ggplot2::aes(label = name, color = labeling, wt = 0.5))
plot1 <- plot1 + ggplot2::ggtitle(paste0("Correlation: ", round(stats::cor(deseq2_lfc, esvd_lfc), 2))) +
  ggplot2::xlab("DESeq2 LFC") + ggplot2::ylab("eSVD-DE LFC")
plot1 <- plot1 + Seurat::NoLegend()

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_esvd-deseq2_lfc.png",
                plot1, device = "png", width = 7, height = 7, units = "in")



