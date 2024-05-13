rm(list=ls())
library(Seurat)

input_folder_path <- "~/kzlinlab/data/sea-ad/"
file_vec <- c("astrocyte_mtg.rds", "microglia_mtg.rds", "oligo_mtg.rds", "opc_mtg.rds")

figure_folder <- "~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5/"
output_folder <- "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/"

# https://satijalab.org/seurat/archive/v4.3/merge
seurat_list <- lapply(file_vec, function(file){
  print(paste0("Working on: ", file))
  seurat_obj <- readRDS(paste0(input_folder_path, file))
  # mat <- Seurat::GetAssayData(object = seurat_obj, 
  #                             assay = "RNA", 
  #                             layer = "data")
  
  return(seurat_obj)
})

seurat_all <- merge(seurat_list[[1]], y = seurat_list[[2]])
seurat_all <- merge(seurat_all, y = seurat_list[[3]])
seurat_all <- merge(seurat_all, y = seurat_list[[4]])

seurat_all
table(seurat_all$Subclass)

rm(ls = seurat_list); gc(TRUE)

print("Finding variable genes")
seurat_all <- Seurat::FindVariableFeatures(seurat_all, 
                                           selection.method = "vst", 
                                           nfeatures = 2000)
seurat_all <- subset(seurat_all, features = Seurat::VariableFeatures(seurat_all))

seurat_all <- Seurat::ScaleData(seurat_all, 
                                features = Seurat::VariableFeatures(seurat_all))

print("Downsampling")
set.seed(10)
n <- length(Seurat::Cells(seurat_all))
keep_idx <- rep(FALSE, n)
keep_idx[sample(1:n, size = round(n/4))] <- TRUE
seurat_all$keep <- keep_idx
seurat_all <- subset(seurat_all, keep == TRUE)

seurat_all

print("Doing PCA")
set.seed(10)
seurat_all <- Seurat::RunPCA(seurat_all, 
                             features = Seurat::VariableFeatures(seurat_all),
                             verbose = FALSE)

print("Doing UMAP")
set.seed(10)
seurat_all <- Seurat::RunUMAP(seurat_all, dims = 1:30)

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(seurat_all, date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_seurat-glial.RData")

print("Done! :)")





