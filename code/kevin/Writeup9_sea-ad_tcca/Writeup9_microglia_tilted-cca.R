# https://linnykos.github.io/tiltedCCA/articles/embryo.html
rm(list=ls())

library(Seurat)
library(tiltedCCA)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_multiome.RData")

set.seed(10)
Seurat::DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                           nfeatures = 2000)
mat_1 <- Matrix::t(SeuratObject::LayerData(seurat_obj,
                                           layer = "data",
                                           assay = "RNA",
                                           features = Seurat::VariableFeatures(seurat_obj)))
Seurat::DefaultAssay(seurat_obj) <- "ATAC"
mat_2 <- Matrix::t(SeuratObject::LayerData(seurat_obj,
                                           layer = "data",
                                           assay = "ATAC"))

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

############

set.seed(10)
multiSVD_obj <- tiltedCCA::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:30, dims_2 = 2:30,
                                           center_1 = T, center_2 = F,
                                           normalize_row = T,
                                           normalize_singular_value = T,
                                           recenter_1 = F, recenter_2 = T,
                                           rescale_1 = F, rescale_2 = T,
                                           scale_1 = T, scale_2 = F,
                                           verbose = 1)

multiSVD_obj <- tiltedCCA::form_metacells(input_obj = multiSVD_obj,
                                          large_clustering_1 = NULL, 
                                          large_clustering_2 = NULL, 
                                          num_metacells = 100,
                                          verbose = 1)

multiSVD_obj <- tiltedCCA::compute_snns(input_obj = multiSVD_obj,
                                        latent_k = 20,
                                        num_neigh = 15,
                                        bool_cosine = T,
                                        bool_intersect = F,
                                        min_deg = 15,
                                        verbose = 1)

multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj,
                                     verbose = 1)

multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                       verbose = 1)

multiSVD_obj <- tiltedCCA::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                   verbose = 1,
                                                   bool_modality_1_full = T,
                                                   bool_modality_2_full = F)

save(multiSVD_obj,
     date_of_run, session_info, 
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_tcca.RData")

set.seed(10)
seurat_obj[["common_tcca"]] <- tiltedCCA::create_SeuratDim(input_obj = multiSVD_obj,
                                                           what = "common",
                                                           aligned_umap_assay = "umap",
                                                           seurat_obj = seurat_obj,
                                                           seurat_assay = "RNA",
                                                           verbose = 1)
set.seed(10)
seurat_obj[["distinct1_tcca"]] <- tiltedCCA::create_SeuratDim(input_obj = multiSVD_obj,
                                                              what = "distinct_1",
                                                              aligned_umap_assay = "umap",
                                                              seurat_obj = seurat_obj,
                                                              seurat_assay = "RNA",
                                                              verbose = 1)
set.seed(10)
seurat_obj[["distinct2_tcca"]] <- tiltedCCA::create_SeuratDim(input_obj = multiSVD_obj,
                                                              what = "distinct_2",
                                                              aligned_umap_assay = "umap",
                                                              seurat_obj = seurat_obj,
                                                              seurat_assay = "RNA",
                                                              verbose = 1)

save(multiSVD_obj, seurat_obj,
     date_of_run, session_info, 
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_tcca.RData")

print("Done! :)")
