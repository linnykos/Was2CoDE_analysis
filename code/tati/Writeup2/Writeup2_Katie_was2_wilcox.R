# Libraries
rm(list=ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_ideascustom.RData")
results_mat <-was2de_pvalue(dist_list, meta_ind)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Wilcoxin_test of the microglia data.",
              "This was done on the data in ~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_ideascustom.RData")
save(results_mat,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_microglia_Katie_was2_wilcox.RData")

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
