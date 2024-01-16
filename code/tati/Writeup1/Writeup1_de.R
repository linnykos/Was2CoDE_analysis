rm(list=ls())

library(Seurat)

load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData")
set.seed(10)

ss_data_norm
head(ss_data_norm@meta.data)
table(ss_data_norm$Pt_ID, ss_data_norm$CognitiveStatus)

Seurat::DefaultAssay(ss_data_norm) <- "integrated"
Seurat::Idents(ss_data_norm) <- "CognitiveStatus"
de_results_wilcoxon <- Seurat::FindMarkers(ss_data_norm, 
                                           features = Seurat::VariableFeatures(ss_data_norm),
                                           test.use = "wilcox",
                                           ident.1 = "Dementia", 
                                           ident.2 = "No_dementia",
                                           assay = "integrated",
                                           slot = "data",
                                           verbose = TRUE)

save(de_results_wilcoxon, 
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup1/Writeup1_de.RData")

de_results_t <- Seurat::FindMarkers(ss_data_norm, 
                                    features = Seurat::VariableFeatures(ss_data_norm),
                                    test.use = "t",
                                    ident.1 = "Dementia", 
                                    ident.2 = "No_dementia",
                                    assay = "integrated",
                                    slot = "data",
                                    verbose = TRUE)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Trying out different DE tests on the data.",
              "Input data from ~/kzlinlab/projects/subject-de/out/kevin/preprocess/processed.RData.")

save(de_results_t,
     de_results_wilcoxon, 
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_de.RData")

print("Done! :)")




