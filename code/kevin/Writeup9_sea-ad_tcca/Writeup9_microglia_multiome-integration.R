rm(list=ls())

library(Seurat)
library(SeuratObject)

seurat_rna <- readRDS("~/kzlinlab/data/sea-ad/microglia_mtg.rds")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_ATACseq_10xMulti_Microglia-PVM.RData")
seurat_atac <- seurat_obj

rna_cells <- Seurat::Cells(seurat_rna)
atac_cells <- Seurat::Cells(seurat_atac)

length(rna_cells)
length(atac_cells)
length(intersect(atac_cells, rna_cells))
cell_names <- intersect(atac_cells, rna_cells)

keep_vec <- rep(FALSE, length(rna_cells))
keep_vec[rna_cells %in% cell_names] <- TRUE
seurat_rna$keep <- keep_vec
seurat_rna <- subset(seurat_rna, keep == TRUE)

keep_vec <- rep(FALSE, length(atac_cells))
keep_vec[atac_cells %in% cell_names] <- TRUE
seurat_atac$keep <- keep_vec
seurat_atac <- subset(seurat_atac, keep == TRUE)

#######

rna_cells <- Seurat::Cells(seurat_rna)
seurat_obj <- seurat_rna

atac_counts <- SeuratObject::LayerData(object = seurat_atac, 
                                       layer = "counts",
                                       assay = "RNA")
atac_counts <- atac_counts[,rna_cells]
seurat_obj[["ATAC"]] <- Seurat::CreateAssayObject(atac_counts)

atac_counts <- SeuratObject::LayerData(object = seurat_atac, 
                                       layer = "data",
                                       assay = "RNA")
atac_counts <- atac_counts[,rna_cells]
SeuratObject::LayerData(object = seurat_obj, 
                        layer = "data",
                        assay = "ATAC") <- atac_counts

######

# port over the embeddings
embedding_mat <- seurat_atac[["python_X_MultiVI"]]@cell.embeddings[rna_cells,]
seurat_obj[["MultiVI"]] <- Seurat::CreateDimReducObject(embedding_mat,
                                                        key = "ATAC")

embedding_mat <- seurat_atac[["python_X_umap"]]@cell.embeddings[rna_cells,]
seurat_obj[["umap_atac"]] <- Seurat::CreateDimReducObject(embedding_mat,
                                                          key = "ATAC")

######

# port over the colorings
seurat_obj@misc <- c(seurat_obj@misc, seurat_atac@misc)

save(seurat_obj,
     date_of_run, session_info, 
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_multiome.RData")


