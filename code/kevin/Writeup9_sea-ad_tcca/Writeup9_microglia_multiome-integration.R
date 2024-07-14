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
seurat_obj[["python_umap_atac"]] <- Seurat::CreateDimReducObject(embedding_mat,
                                                          key = "ATAC")

seurat_obj[["python_umap_rna"]] <- seurat_obj[["umap"]]
seurat_obj[["umap"]] <- NULL

######

# port over the colorings
seurat_obj@misc <- c(seurat_obj@misc, seurat_atac@misc)

###########################

# fix covariates
variable <- "Braak stage"
tmp <- seurat_obj@meta.data[,variable]
tmp2 <- sapply(as.character(tmp), function(x){
  strsplit(x, split = " ")[[1]][2]
})
tmp2[tmp2 == "I"] <- "1"; tmp2[tmp2 == "II"] <- "2"
tmp2[tmp2 == "III"] <- "3"; tmp2[tmp2 == "IV"] <- "4"
tmp2[tmp2 == "V"] <- "5"; tmp2[tmp2 == "VI"] <- "6"
seurat_obj@meta.data[,variable] <- as.numeric(tmp2)

variable <- "Thal phase"
tmp <- seurat_obj@meta.data[,variable]
tmp2 <- sapply(as.character(tmp), function(x){
  strsplit(x, split = " ")[[1]][2]
})
seurat_obj@meta.data[,variable] <- as.numeric(tmp2)

variable <- "LATE-NC stage"
tmp <- seurat_obj@meta.data[,variable]
levels(tmp)[which(levels(tmp) == "Staging Precluded by FTLD with TDP43 or ALS/MND or TDP-43 pathology is unclassifiable")] <- "Unclassifiable"
seurat_obj@meta.data[,variable] <- tmp

variable <- "Age at death"
tmp <- seurat_obj@meta.data[,"development_stage"]
tmp <- as.character(tmp)
tmp[which(tmp == "80 year-old and over human stage")] <- "90-year-old human stage"
tmp2 <- sapply(tmp, function(x){
  strsplit(x, split = "-")[[1]][1]
})
seurat_obj@meta.data[,variable] <- as.numeric(tmp2)



#################

# massive overhaul to remove spaces

df <- seurat_obj@meta.data
colnames(df) <- sapply(colnames(df), function(x){
  gsub(pattern = " ", replacement = "_", x = x)
})

for(j in 1:ncol(df)){
  if(is.factor(df[,j]) | is.character(df[,j])){
    print(paste("Processing", colnames(df)[j]))
    
    vec <- df[,j]
    factor_bool <- FALSE
    if(is.factor(vec)) {
      factor_bool <- TRUE
      vec <- as.character(vec)
    }
    vec <- gsub(pattern = " ", replacement = "_", x = vec)
    if(factor_bool) vec <- factor(vec)
    df[,j] <- vec
  }
}

seurat_obj@meta.data <- df

summary(seurat_obj@meta.data)

#################

# create color gradient for donors
tab_mat <- table(seurat_obj$donor_id, seurat_obj@meta.data[,"Cognitive_status"])
tab_mat <- tab_mat[,c("No_dementia", "Dementia")]
tab_mat <- tab_mat[rowSums(tab_mat)!=0,]

ad_donor <- rownames(tab_mat)[which(tab_mat[,"Dementia"] != 0)]
ctrl_donor <- rownames(tab_mat)[which(tab_mat[,"Dementia"] == 0)]

case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(length(ad_donor))
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(length(ctrl_donor))
col_vec <- c(case_color_palette, control_color_palette)
names(col_vec) <- c(ad_donor, ctrl_donor)

seurat_obj@misc$donor_id_colors <- col_vec

##########################

# make some clusterings

set.seed(10)
Seurat::DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, 
                                           selection.method = "vst", nfeatures = 2000)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             features = Seurat::VariableFeatures(seurat_obj),
                             verbose = FALSE)
set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)

set.seed(10)
Seurat::DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- Signac::RunTFIDF(seurat_obj)
seurat_obj <- Signac::FindTopFeatures(seurat_obj, min.cutoff = "q0")
seurat_obj <- Signac::RunSVD(seurat_obj)
set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    reduction = 'lsi',
                                    dims = 2:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "lsi", 
                              dims = 2:30, 
                              reduction.name = "umap_atac", 
                              reduction.key = "atacUMAP_")

set.seed(10)
seurat_obj <- Seurat::FindMultiModalNeighbors(
  seurat_obj, 
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:30), 
  modality.weight.name = "RNA.weight"
)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, nn.name = "weighted.nn", 
                              reduction.name = "umap_wnn", 
                              reduction.key = "wnnUMAP_")

###########################

save(seurat_obj,
     date_of_run, session_info, 
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_multiome.RData")


