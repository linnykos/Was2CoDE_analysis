rm(list=ls())
library(openxlsx)
library(tibble)
library(dplyr)

metadata <- openxlsx::read.xlsx("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_cohort_donor_metadata_082222.xlsx")
neuropath <- read.csv("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv")
neuropath <- neuropath[,-1]

all(metadata$Donor.ID %in% neuropath$Donor.ID)

metadata2 <- tibble::as_tibble(metadata)
neuropath2 <- tibble::as_tibble(neuropath)

tib <- dplyr::left_join(
  x = metadata2,
  y = neuropath2,
  by = "Donor.ID"
)

colname_vec <- colnames(tib)

#######################

# APOE4.Status, Sex
# Age.at.Death, Overall.AD.neuropathological.Change, Thal, Braak, CERAD.score, Overall.CAA.Score, LATE, Total.microinfarcts.in.screening.sections






