rm(list=ls())
library(EnhancedVolcano)
library(openxlsx)
library(ggplot2)
library(IdeasCustom)
library(dplyr)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_NEBULA.RData")

head(nebula_res$summary)
res <- nebula_res$summary

pval_vec <- res[,"p_Study_DesignationAD"]
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)
#define thresholds
pCutoff <- max(pval_vec[idx])
FCcutoff <- quantile(abs(res[,"logFC_Study_DesignationAD"]), probs = 0.9)
xlim <- c(-1,1) * quantile(abs(res[,"logFC_Study_DesignationAD"]), probs = 0.99)

idx <- which(abs(res[,"logFC_Study_DesignationAD"]) <= max(xlim))
ylim <- c(0, 20)

plot1 <- EnhancedVolcano::EnhancedVolcano(
  res,
  lab = res$gene,
  x = "logFC_Study_DesignationAD",
  y = "p_Study_DesignationAD",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  xlim = xlim,
  ylim = ylim
)

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/Writeup5_prater_NEBULA_volcano.png",
                plot1, device = "png", width = 7, height = 7, units = "in")

#########################
res <- nebula_res$summary

microglia_prater <- unique(microglia_prater_df$Gene[microglia_prater_df$DEG])
housekeeping_hounkpe <- housekeeping_hounkpe_df$Gene

head(nebula_res$summary)
res <- nebula_res$summary
#define thresholds
pval_vec <- res[,"p_Study_DesignationAD"]
pval_adj_vec <- stats::p.adjust(pval_vec, method = "BH")
idx <- which(pval_adj_vec <= 0.05)

pCutoff <- max(pval_vec[idx])
FCcutoff <- quantile(abs(res[,"logFC_Study_DesignationAD"]), probs = 0.9)
xlim <- c(-1,1) * quantile(abs(res[,"logFC_Study_DesignationAD"]), probs = 0.99)
ylim <- c(0, 20)

res$GeneType <- ifelse(res$gene %in% microglia_prater, "Prater",
                       ifelse(res$gene %in% housekeeping_hounkpe, "Housekeeping", "Other"))

# First, reorder the levels of GeneType
res$GeneType <- factor(res$GeneType, levels = c("Other", "Housekeeping", "Prater"))
res <- res[c(which(res$GeneType == "Other"),
             which(res$GeneType == "Housekeeping"),
             which(res$GeneType == "Prater")),
           ]

# Now create the plot
plot0 <- ggplot(res, aes(x = logFC_Study_DesignationAD, y = -log10(p_Study_DesignationAD), color = GeneType)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Prater" = "red", "Housekeeping" = "green", "Other" = "grey")) +
  geom_text(aes(label = ifelse(gene %in% c(microglia_prater, housekeeping_hounkpe), as.character(gene), "")),
            vjust = 1.5, hjust = 0.5, check_overlap = TRUE, size = 3) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme_minimal() +
  xlim(xlim) + ylim(ylim) + 
  geom_vline(xintercept = c(-FCcutoff, FCcutoff), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -log10(pCutoff), linetype = "dashed", color = "black")

ggplot2::ggsave(filename = "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/Writeup5_prater_NEBULA_volcano_colored_by_GeneLists.png",
                plot0, device = "png", width = 7, height = 7, units = "in")

genes_above_threshold <- res %>%
  filter(-log10(p_Study_DesignationAD) > -log10(pCutoff)) %>%
  select(gene)

print(genes_above_threshold)

write.csv(genes_above_threshold, "~/kzlinlab/projects/subject-de/out/tati/Writeup5/genes_above_threshold_NEBULA.csv", row.names = FALSE)
write.table(genes_above_threshold, "~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup5/genes_above_threshold_NEBULA.csv", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
#########################
# sun_sheet <- openxlsx::read.xlsx(
#   xlsxFile = "../../../../../data/1-s2.0-S0092867423009716-mmc1.xlsx",
#   sheet = "Page 10.DEGs_AD"
# )
# sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
# bool_vec <- rep(FALSE, nrow(res))
# names(bool_vec) <- res$gene
# bool_vec[which(res$gene %in% sun_genes)] <- TRUE
# table(bool_vec)
# 
# # create a custom volcano plot
# 
# lfc_vec <- res[,"logFC_Study_DesignationAD"]
# xlim <- c(-1,1) * quantile(abs(lfc_vec), probs = 0.99)
# 
# pvalue_vec <- res[,"p_Study_DesignationAD"]
# logpval_vec <- -log10(pvalue_vec)
# df <- data.frame(lfc = lfc_vec,
#                  log10pval = logpval_vec,
#                  name = res$gene,
#                  labeling = bool_vec)
# # put all the labeling == TRUE on bottom
# df <- df[c(which(df[,"labeling"] == "FALSE"), which(df[,"labeling"] == "TRUE")),]
# 
# plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lfc,
#                                           y = log10pval))
# plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
# plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
# plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == "TRUE"),
#                                           ggplot2::aes(label = name, color = labeling),
#                                           box.padding = ggplot2::unit(0.5, 'lines'),
#                                           point.padding = ggplot2::unit(1.6, 'lines'),
#                                           max.overlaps = 50)
# plot1 <- plot1 + ggplot2::xlim(xlim)
# if(any(pvalue_vec <= 0.05)) {
#   plot1 <- plot1 + ggplot2::geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed",
#                                        color = "red", linewidth=2)
# }
# plot1 <- plot1 + ggplot2::ggtitle(paste0("NEBULA volcano plot")) +
#   ggplot2::xlab("Log fold change") + ggplot2::ylab("NEBULA p-value (-Log10)")
# plot1 <- plot1 + Seurat::NoLegend()
# 
# ggplot2::ggsave(filename = "../../../figures/kevin/Writeup1/Writeup1_nebula_volcano-annotated.png",
#                 plot1, device = "png", width = 5, height = 5, units = "in")
