rm(list=ls())
library(Seurat)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")

umap_list <- sapply(1:10, function(i){
  print(paste("Working on trial:", i))
  seurat_all[["umap"]] <- NULL
  
  set.seed(i)
  # see https://jlmelville.github.io/uwot/articles/abparams.html
  seurat_all <- Seurat::RunUMAP(seurat_all, 
                                dims = 1:30,
                                seed.use = i,
                                min.dist = 0.1, #default: 0.3
                                spread = 0.1) #default: 1
  return(seurat_all[["umap"]])
})

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(umap_list, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_glial_umaps_alt2.RData")

print("Done! :)")