# https://github.com/mojaveazure/seurat-disk/issues/166
# We need the older version of Seurat
# remotes::install_version("SeuratObject", 
#                          version = "4.1.4", 
#                          repos = c("https://satijalab.r-universe.dev", getOption("repos")),
#                          lib = "/home/users/kzlin/R/x86_64-pc-linux-gnu-library/Seuratv4")
# remotes::install_version("Seurat", 
#                          version = "4.4.0", 
#                          repos = c("https://satijalab.r-universe.dev", getOption("repos")),
#                          lib = "/home/users/kzlin/R/x86_64-pc-linux-gnu-library/Seuratv4")
# remotes::install_github("mojaveazure/seurat-disk",
#                         force = T,
#                         lib = "/home/users/kzlin/R/x86_64-pc-linux-gnu-library/Seuratv4")

rm(list=ls())
library(SeuratObject, lib.loc = "/home/users/kzlin/R/x86_64-pc-linux-gnu-library/Seuratv4")
library(Seurat, lib.loc = "/home/users/kzlin/R/x86_64-pc-linux-gnu-library/Seuratv4")
library(SeuratDisk, lib.loc = "/home/users/kzlin/R/x86_64-pc-linux-gnu-library/Seuratv4")

set.seed(10)
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess_Seuratv3.RData")

# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
# https://github.com/mojaveazure/seurat-disk/issues/166
SeuratDisk::SaveH5Seurat(seurat_diet, 
                         filename = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia.h5Seurat",
                         overwrite = T)
SeuratDisk::Convert("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia.h5Seurat", 
                    dest = "h5ad",
                    overwrite = TRUE)