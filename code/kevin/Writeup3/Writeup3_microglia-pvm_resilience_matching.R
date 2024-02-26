rm(list=ls())
library(openxlsx)
library(optmatch)
source("matching.R")

metadata <- openxlsx::read.xlsx("~/kzlinlab/data/sea-ad/sea-ad_cohort_donor_metadata_082222.xlsx")
rownames(metadata) <- metadata$Donor.ID
harmonized_scores <- openxlsx::read.xlsx("~/kzlinlab/data/sea-ad/joey-harmonized-cognitive/sea_ad_cognitive_slopes_79.xlsx")
rownames(harmonized_scores) <- harmonized_scores$donor_name
neuropath <- read.csv("~/kzlinlab/data/sea-ad/sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv")
rownames(neuropath) <- neuropath$Donor.ID
path_idx <- grep("^percent.*area_Grey", colnames(neuropath))
colnames(neuropath)[path_idx]
neuropath <- neuropath[,path_idx]

included_idx <- which(metadata$Donor.ID %in% harmonized_scores$donor_name)
metadata <- metadata[included_idx,]

ad_idx <- intersect(
  which(metadata$Braak %in% c("Braak V", "Braak VI")),
  which(metadata$CERAD.score %in% c("Moderate", "Frequent"))
)
metadata <- metadata[ad_idx,]
harmonized_scores <- harmonized_scores[rownames(metadata),]
neuropath <- neuropath[rownames(metadata),]

col_vec <- colnames(neuropath)
col_vec <- sapply(col_vec, function(x){
  strsplit(x, split = "\\.")[[1]][2]
})
colnames(neuropath) <- col_vec
neuropath <- scale(neuropath)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
esvd_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat) + tcrossprod(eSVD_obj$covariates[,"resilient_Resilient"], eSVD_obj$fit_Second$z_mat[,"resilient_Resilient"])

keep_vec <- rep(FALSE, length(SeuratObject::Cells(seurat_obj)))
keep_vec[which(SeuratObject::Cells(seurat_obj) %in% rownames(esvd_mat))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)
cell_vec <- SeuratObject::Cells(seurat_obj)

# process the seurat_obj$donor_id
donor_vec <- seurat_obj$donor_id
donor_vec <- sapply(donor_vec, function(x){
  paste0(substr(x, start = 0, stop = 3), ".", 
         substr(x, start = 4, stop = 5), ".",
         substr(x, start = 6, stop = 8))
})
seurat_obj$donor_id <- factor(donor_vec)

donor_list <- lapply(rownames(neuropath), function(donor){
  tmp <- which(seurat_obj$donor_id == donor)
  names(tmp) <- NULL
  tmp
})
names(donor_list) <- rownames(neuropath)
donor_numcell <- sapply(donor_list, length)

if(any(donor_numcell == 0)){
  keep_donor <- names(donor_list)[which(donor_numcell > 0)]
  
  metadata <- metadata[keep_donor,]
  harmonized_scores <- harmonized_scores[keep_donor,]
  neuropath <- neuropath[keep_donor,]
  donor_list <- donor_list[keep_donor]
}

seurat_obj$donor_id <- factor(donor_vec)

#################

cognition_vec <- harmonized_scores[,"slope_zmem0"]
med_val <- stats::median(cognition_vec)
cognition_vec2 <- rep(FALSE, length(cognition_vec))
cognition_vec2[cognition_vec > med_val] <- TRUE
df <- data.frame(neuropath)
df$cognition <- cognition_vec2

matching_res <- matching_func(
  df = df,
  class_covariate = "cognition"
)

######################

# for a given gene:
# gene <- "SLC2A3"
gene_names <- colnames(esvd_mat)

avg_mat <- sapply(donor_list, function(idx_vec){
  Matrix::colMeans(esvd_mat[idx_vec,])
})
avg_mat <- t(avg_mat)

res_mat <- sapply(gene_names, function(gene){
  compute_wilcoxon(bg_matches = matching_res$bg_matches,
                   gene = gene,
                   avg_mat = avg_mat,
                   signal_matches = matching_res$signal_matches)
})

res_df <- data.frame(
  pval = as.numeric(res_mat["pval",]),
  side = res_mat["side",]
)
rownames(res_df) <- colnames(res_mat)

# sun_sheet <- openxlsx::read.xlsx(
#   xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
#   sheet = "Page 10.DEGs_AD"
# ) 
# sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
# sun_genes <- intersect(gene_names, sun_genes)
# 
# hk_df <- read.csv("~/kzlinlab/projects/subject-de/data/Housekeeping_GenesHuman.csv", sep = ";")
# hk_genes <- hk_df[,"Gene.name"]
# hk_genes <- intersect(gene_names, hk_genes)

# res_mat2[,which(colnames(res_mat2) %in% sun_genes)]
# pvalue_vec <- res_mat["twosided",]
# candidate_genes <- colnames(res_mat)[order(abs(pvalue_vec - 0.001), decreasing = F)[1:100]]
# candidate_genes[candidate_genes %in% sun_genes]
# candidate_genes[candidate_genes %in% hk_genes]

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(res_df, matching_res,
     df, avg_mat, 
     date_of_run, session_info, 
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching.RData")

