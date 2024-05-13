rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")

umap_list <- sapply(1:10, function(i){
  print(paste("Working on trial:", i))
  seurat_all[["umap"]] <- NULL
  
  set.seed(i)
  seurat_all <- Seurat::RunUMAP(seurat_all, 
                                dims = 1:30,
                                seed.use = i)
  return(seurat_all[["umap"]])
})

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(umap_list, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_glial_umaps.RData")

print("Done! :)")