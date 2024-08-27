rm(list=ls())

library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData")
set.seed(10)

ss_data_norm
head(ss_data_norm@meta.data)
table(ss_data_norm$Pt_ID, ss_data_norm$CognitiveStatus)


# INSERT CODE TO RUN IDEAS AND GET THE P VALUES

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save( # SAVE THE RESULTS 
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_de.RData")

print("Done! :)")




