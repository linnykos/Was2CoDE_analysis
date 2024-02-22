rm(list=ls())
library(Seurat)
library(ideas)
library(foreach)
library(doRNG)
RNGkind("L'Ecuyer-CMRG")
library(doParallel)
doParallel::registerDoParallel(cores=6)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_ideas.RData")
set.seed(10)

neuropath <- read.csv("~/kzlinlab/data/sea-ad/sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv")
rownames(neuropath) <- neuropath$Donor.ID
path_idx <- grep("^percent.*area_Grey", colnames(neuropath))
neuropath <- neuropath[,path_idx]
neuropath <- neuropath[indiv_covariates$individual,]
neuropath <- scale(neuropath)
indiv_covariates <- cbind(indiv_covariates, neuropath)

print("Starting ideas::ideas_dist")
set.seed(10)
time_start <- Sys.time()
pval_res <- ideas::permanova(dist_tensor,
                             meta_ind = indiv_covariates,
                             var2test = "resiliency",
                             var2adjust = colnames(),
                             var2test_type = "binary",
                             n_perm = 4999,
                             r.seed = 904)
time_end <- Sys.time()

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(date_of_run, session_info, 
     time_start, time_end,
     cell_covariates, dist_tensor, indiv_covariates,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_ideas.RData")

print("Done! :)")