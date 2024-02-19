rm(list=ls())
library(openxlsx)
library(tibble)
library(dplyr)

metadata <- openxlsx::read.xlsx("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_cohort_donor_metadata_082222.xlsx")
neuropath <- read.csv("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv")
neuropath <- neuropath[,-1]

all(metadata$Donor.ID %in% neuropath$Donor.ID)
# 
# metadata2 <- tibble::as_tibble(metadata)
# neuropath2 <- tibble::as_tibble(neuropath)
# 
# tib <- dplyr::left_join(
#   x = metadata2,
#   y = neuropath2,
#   by = "Donor.ID"
# )

colname_vec <- colnames(neuropath)

perctange_idx <- grep("perc", colname_vec)
colname_vec[perctange_idx]

#######################

# preprocess certain columns
subj_idx <- intersect(
  intersect(which(metadata$Braak %in% c("Braak V", "Braak VI")),
            which(metadata$CERAD.score %in% c("Moderate", "Frequent"))),
  which(metadata$Interval.from.last.CASI.in.months <= 36)
)
subj_id <- metadata$Donor.ID[subj_idx]

for(i in perctange_idx){
  hist(neuropath[which(neuropath$Donor.ID %in% subj_id),i], 
       col = "gray",
       main = colnames(neuropath)[i])
}


# APOE4.Status, Sex
# Age.at.Death, Overall.AD.neuropathological.Change, Thal, Braak, CERAD.score, Overall.CAA.Score, LATE, Total.microinfarcts.in.screening.sections






