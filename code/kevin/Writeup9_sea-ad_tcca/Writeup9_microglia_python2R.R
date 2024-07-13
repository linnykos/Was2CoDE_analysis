rm(list=ls()); gc(TRUE)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
options(Seurat.object.assay.version = "v5")

# from https://github.com/theislab/zellkonverter/issues/38
adata <- zellkonverter::readH5AD(
  file = "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_SEAAD_MTG_ATACseq_10xMulti_Microglia-PVM.h5ad"
)

names(adata@assays)

tmp <- Seurat::as.Seurat(adata, counts = "X", data = "X")
# converted most of the things. It creates an assay by default called "originalexp"
# see also https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html

seurat_obj <- Seurat::CreateSeuratObject(
  counts = tmp[["originalexp"]],
  data = tmp[["originalexp"]],
  meta.data = tmp@meta.data
)

seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

# put in the assays
name_vec <- names(adata@assays)
gene_vec <- SeuratObject::Features(seurat_obj[["RNA"]])

# put in the gene metafeatures
gene_metadata <- SingleCellExperiment::rowData(adata)
seurat_obj[["RNA"]]@misc <- as.data.frame(gene_metadata)

# put in the dimension reductions
name_vec <- SingleCellExperiment::reducedDimNames(adata)
for(name_val in name_vec){
  mat <- SingleCellExperiment::reducedDim(adata, name_val)
  name_val2 <- paste0("python_", name_val)
  colnames(mat) <- paste0(name_val2, "_", 1:ncol(mat))
  
  seurat_obj[[name_val2]] <- Seurat::CreateDimReducObject(embeddings = mat,
                                                          assay = "RNA")
}

# put in the metadata
metadata_list <- adata@metadata
idx <- which(sapply(1:length(metadata_list), function(i){
  class(metadata_list[[i]]) %in% c("dgCMatrix", "dgRMatrix")
}))
if(length(idx) > 0){
  graph_list <- metadata_list[idx]
  metadata_list <- metadata_list[-idx]
} else {
  graph_list <- numeric(0)
}
seurat_obj@misc <- metadata_list

if(length(graph_list) > 0){
  for(name_val in names(graph_list)){
    print(paste0("Putting in graph ", name_val))
    
    seurat_obj@graphs[[name_val]] <- graph_list[[name_val]]
    rownames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)
    colnames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)
    
  }
}

Seurat::DefaultAssay(seurat_obj) <- "RNA"

# for the cluster color specifically, add the names for convenience
names(seurat_obj@misc[["APOE4.Status_colors"]]) <- sort(unique(seurat_obj$APOE.Genotype))
names(seurat_obj@misc[["Braak_colors"]]) <- sort(unique(seurat_obj$Braak))
names(seurat_obj@misc[["CERAD.score_colors"]]) <- sort(unique(seurat_obj$CERAD.score))
names(seurat_obj@misc[["Cognitive.Status_colors"]]) <- sort(unique(seurat_obj$Cognitive.Status))
names(seurat_obj@misc[["Highest.Lewy.Body.Disease_colors"]]) <- sort(unique(seurat_obj$Highest.Lewy.Body.Disease))
names(seurat_obj@misc[["LATE_colors"]]) <- sort(unique(seurat_obj$LATE))
names(seurat_obj@misc[["Overall.AD.neuropathological.Change_colors"]]) <- sort(unique(seurat_obj$Overall.AD.neuropathological.Change))
names(seurat_obj@misc[["Sex_colors"]]) <- sort(unique(seurat_obj$Sex))
names(seurat_obj@misc[["Subclass_colors"]]) <- sort(unique(seurat_obj$Subclass))
names(seurat_obj@misc[["Thal_colors"]]) <- sort(unique(seurat_obj$Thal))

# clearing environment
ls_vec <- ls()
ls_vec <- ls_vec[ls_vec != "seurat_obj"]
rm(list = ls_vec)
gc(TRUE)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Conversion of SEA-AD microglia PVM's ATAC data into R")

save(seurat_obj,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_ATACseq_10xMulti_Microglia-PVM.RData")
