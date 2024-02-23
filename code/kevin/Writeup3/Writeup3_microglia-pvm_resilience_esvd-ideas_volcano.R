rm(list=ls())
library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_ideas_pvalues.RData")

set.seed(10)

mat <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat) + tcrossprod(eSVD_obj$covariates[,"resilient_Resilient"], eSVD_obj$fit_Second$z_mat[,"resilient_Resilient"])

resiliency_vec <- which(eSVD_obj$covariates[,"resilient_Resilient"] == 1)
nonresilient_mean <- Matrix::colMeans(mat[-resiliency_vec,])
resilient_mean <- Matrix::colMeans(mat[resiliency_vec,])
lfc_vec <- resilient_mean - nonresilient_mean

pvalue_vec <- pval_res[names(lfc_vec)]

pval_adj_vec <- stats::p.adjust(pvalue_vec, method = "BH")
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

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_ideas_volcano.png",
                plot1, device = "png", width = 7, height = 7, units = "in")


#########################

logpval_vec <- -log10(pvalue_vec)
fdr_vec <- pval_adj_vec
fdr_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
fdr_bool <- rep(FALSE, length(fdr_vec))
names(fdr_bool) <- names(fdr_vec)
fdr_bool[which(fdr_vec <= 0.05)] <- TRUE

custom_volcano_highlight_func <- function(
    custom_genes,
    fdr_bool,
    lfc_vec, 
    logpval_vec,
    main_addition
){
  gene_vec <- names(lfc_vec)
  fdr_genes <- gene_vec[which(fdr_bool == TRUE)]
  
  bool_vec <- rep(FALSE, length(fdr_bool))
  names(bool_vec) <- names(fdr_bool)
  bool_vec[which(names(bool_vec) %in% custom_genes)] <- TRUE
  bool_vec1 <- bool_vec*fdr_bool # both in validation dataset and our FDR
  bool_vec2 <- bool_vec - bool_vec1 # only in validation
  label_vec <- rep("None", length(bool_vec))
  label_vec[which(bool_vec2 == 1)] <- "Validation"
  label_vec[which(bool_vec1 == 1)] <- "Both"
  
  m <- length(custom_genes)
  n <- length(gene_vec) - m
  k <- length(fdr_genes)
  x <- length(intersect(fdr_genes, custom_genes))
  fisher <- sum(sapply(x:k, function(i){
    stats::dhyper(x = i, m = m, n = n, k = k, log = F)
  }))
  
  # create a custom volcano plot 
  df <- data.frame(lfc = lfc_vec,
                   log10pval = logpval_vec,
                   name = gene_vec,
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
  plot1 <- plot1 + ggplot2::ggtitle(paste0("IDEAS volcano plot, ", main_addition,
                                           "\nFisher -log10pvalue: ", round(-log10(fisher), 2))) +
    ggplot2::xlab("Log fold change") + ggplot2::ylab("IDEAS p-value (-Log10)")
  plot1 <- plot1 + Seurat::NoLegend()
  
  plot1
}

#########################

prater_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/43587_2023_424_MOESM3_ESM.xlsx",
  sheet = "AD_vs_Ctrl_DEGs",
  startRow = 3
) 
prater_genes <- sort(unique(prater_sheet[which(prater_sheet[,"padj"] <= 0.05),"Gene"]))

plot1 <- custom_volcano_highlight_func(
  custom_genes = prater_genes,
  fdr_bool = fdr_bool,
  lfc_vec = lfc_vec, 
  logpval_vec = logpval_vec,
  main_addition = "Prater genes"
)
ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_ideas_volcano_annotated-prater.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

#########################

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))

plot1 <- custom_volcano_highlight_func(
  custom_genes = sun_genes,
  fdr_bool = fdr_bool,
  lfc_vec = lfc_vec, 
  logpval_vec = logpval_vec,
  main_addition = "Sun genes"
)
ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_ideas_volcano_annotated-sun.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

#########################

mostafavi_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/41593_2018_154_MOESM4_ESM.xlsx",
  sheet = "TableS2_gene_trait_associations",
  startRow = 3
) 
tmp <- 10^(-abs(mostafavi_sheet[,"ClinAD"]))
tmp <- stats::p.adjust(tmp, method = "BH")
mostafavi_genes <- mostafavi_sheet[which(tmp <= 0.05),"Gene.Symbol"]

plot1 <- custom_volcano_highlight_func(
  custom_genes = mostafavi_genes,
  fdr_bool = fdr_bool,
  lfc_vec = lfc_vec, 
  logpval_vec = logpval_vec,
  main_addition = "Mostafavi genes"
)
ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_ideas_volcano_annotated-mostafavi.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###############

brase_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/41467_2023_37437_MOESM13_ESM.xlsx",
  sheet = "h. Microglia_ADADvsCO"
) 
brase_genes <- brase_sheet[,1]

plot1 <- custom_volcano_highlight_func(
  custom_genes = brase_genes,
  fdr_bool = fdr_bool,
  lfc_vec = lfc_vec, 
  logpval_vec = logpval_vec,
  main_addition = "Brase genes"
)
ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_ideas_volcano_annotated-brase.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

############################

# plot the cross-method p-value

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_ideas_pvalues.RData")
gene_vec <- names(pval_res)
ideas_pvalue <- pval_res
ideas_log10 <- -log10(ideas_pvalue)
ideas_fdr <- stats::p.adjust(ideas_pvalue, method = "BH")
ideas_cutoff <- min(ideas_log10[which(ideas_fdr <= 0.05)])

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
label_vec[intersect(which(ideas_log10 > ideas_cutoff),
                    which(gene_vec %in% sun_genes))] <- "Both"
label_vec[intersect(which(esvd_log10 > esvd_cutoff),
                    which(gene_vec %in% sun_genes))] <- "Both"

# create a custom volcano plot 
df <- data.frame(ideas_log10 = ideas_log10,
                 esvd_log10 = esvd_log10,
                 name = gene_vec,
                 labeling = factor(label_vec))
# put all the labeling == TRUE on bottom
df <- df[c(which(df[,"labeling"] == "None"), 
           which(df[,"labeling"] == "Validation"),
           which(df[,"labeling"] == "Both")),]

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = ideas_log10, 
                                          y = esvd_log10))
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c(Both = "red",
                                                       None = "black", 
                                                       Validation = "cadetblue2"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "Both"),
                                          ggplot2::aes(label = name, color = labeling, wt = 0.5))
plot1 <- plot1 + ggplot2::geom_hline(yintercept=esvd_cutoff, linetype="dashed", 
                                     color = "red", linewidth=2)
plot1 <- plot1 + ggplot2::geom_vline(xintercept=ideas_cutoff, linetype="dashed", 
                                     color = "red", linewidth=2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Correlation: ", round(stats::cor(ideas_log10, esvd_log10), 2))) +
  ggplot2::xlab("IDEAS p-value (-Log10)") + ggplot2::ylab("eSVD-DE p-value (-Log10)")
plot1 <- plot1 + Seurat::NoLegend()

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_esvd-ideas_volcano.png",
                plot1, device = "png", width = 7, height = 10, units = "in")


