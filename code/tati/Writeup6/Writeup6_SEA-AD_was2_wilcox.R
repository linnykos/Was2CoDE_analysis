# Libraries
rm(list=ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)
library(IdeasCustom)
set.seed(10)

load("~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_microglia_ideascustom.RData")
# sum(is.na(meta_ind))
meta_ind$ADNC <- with(meta_ind, 
                        ifelse(ADNC %in% c("NotAD", "Low"), "Control", 
                               ifelse(ADNC %in% c("Intermediate", "High"), "Case", ADNC)))

meta_ind$Dementia_Status <- ifelse(meta_ind$ADNC == "Control", "No_dementia", "Dementia")
# table(meta_ind$Dementia_Status)
results_mat <- IdeasCustom::was2de_pvalue(dist_list, meta_ind, "Dementia_Status", verbose = 1)
# summary(results_mat)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Wilcoxin_test of the microglia data.",
              "This was done on the data in ~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_microglia_ideascustom.RData")

save(results_mat,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_was2_wilcox.RData")

print("Done! :)")

# 
# x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
# y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
# 
# mean(x); mean(y)
# 
# stats::wilcox.test(x, y)
# 
# stats::wilcox.test(x, y, alternative = "less") 
# # HA: x < y (we would expect a large p-value, i.e., not significant)
# 
# 
# stats::wilcox.test(x, y, alternative = "greater")  
# # HA: x > y (we would expect a small p-value, i.e., significant)
