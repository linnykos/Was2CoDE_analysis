rm(list=ls())

library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)

df <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/out/tati/Writeup5/Prater_dataset_ingredients.csv")
rownames(df) <- df$Gene
df <- df[order(df$Was2_pval, decreasing = FALSE),]

was2_pval <- df[,"Was2_pval"]
was2_pval <- -log10(was2_pval)
names(was2_pval) <- rownames(df)

teststat_vec <- sort(was2_pval, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500,            # maximum gene set size
  scoreType = "pos"
)

gse_df <- as.data.frame(gse)

gse_df[,c("Description", "NES")]

###########


Was2_logFC <- df[,"Was2_logFC"]
names(Was2_logFC) <- rownames(df)

teststat_vec <- sort(Was2_logFC, decreasing = TRUE)

set.seed(10)
gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)

gse_df <- as.data.frame(gse)

gse_df[,c("Description", "NES")]

pval_adj <- stats::p.adjust(df[,"Was2_pval"], method = "BH")
idx <- which(pval_adj <= 0.05)
fdr_genes <- names(was2_pval)[idx]
color_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(pval_adj))
names(color_vec) <- names(was2_pval)
color_vec[fdr_genes] <- "red"

plot(x = Was2_logFC,
     y = was2_pval,
     col = color_vec,
     pch = 16)

sort(Was2_logFC, decreasing = TRUE)[1:5]
names(sort(Was2_logFC, decreasing = TRUE)[1:5])

##########

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/out/tati/Writeup5/Writeup5_microglia_ideascustom.RData")

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

donor_df <- ss_data_norm@meta.data[,c("Pt_ID","Study_Designation")]
donor_df <- unique(donor_df)
case_pt <- donor_df[which(donor_df$Study_Designation == "AD"),"Pt_ID"]
control_pt <- donor_df[which(donor_df$Study_Designation == "Ctrl"),"Pt_ID"]

case_pt <- as.character(case_pt)
control_pt <- as.character(control_pt)

# https://github.com/TatiZhang/IdeasCustom/blob/main/R/divergence.R
.compute_mean_difference <- function(mat,
                                     case_pt,
                                     control_pt){
  stopifnot(all(rownames(mat) == colnames(mat)),
            length(rownames(mat)) == length(colnames(mat)))
  
  diag(mat) <- NA
  case_idx <- which(rownames(mat) %in% case_pt)
  control_idx <- which(rownames(mat) %in% control_pt)
  
  case_mean <- mean(mat[case_idx,case_idx], na.rm = TRUE)
  control_mean <- mean(mat[control_idx,control_idx], na.rm = TRUE)
  
  case_mean - control_mean
}

logfc_vec <- sapply(1:dim(dist_list[["location_sign"]])[1], function(j){
  mat <- dist_list[["location_sign"]][j,,]
  .compute_mean_difference(
    mat = mat,
    case_pt = case_pt,
    control_pt = control_pt
  )
})

plot(x = logfc_vec,
     y = was2_pval)

# we probably need to do a signed thing -- so plot the % mean difference, and then sign it?


