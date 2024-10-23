rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(arrow)

out_folder <- "~/kzlinlab/projects/subject-de/out/kevin/Writeup11/"
load(paste0(out_folder, "Writeup11_rosmap_scVI-postprocessed.RData"))

metadata <- seurat_obj@meta.data
colnames(metadata)
metadata <- metadata[which(metadata$ADpath == "yes"),]
metadata <- metadata[which(!is.na(metadata$cts_mmse30_lv)),]
metadata <- metadata[,c("Pt_ID", "age_death", "sex")]
median(table(metadata$Pt_ID))
metadata <- unique(metadata)
summary(metadata)

############

rm(list=ls())
library(openxlsx)
library(optmatch)

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

metadata <- metadata[,c("Donor.ID", "Age.at.Death", "Sex")]
tmp <- metadata$Age.at.Death 
tmp[tmp == "90+"] <- "90"
metadata$Age.at.Death <- as.numeric(tmp)
metadata$Sex <- factor(metadata$Sex)
summary(metadata)



