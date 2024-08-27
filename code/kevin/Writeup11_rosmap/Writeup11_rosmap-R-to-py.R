rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)

out_folder <- "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup11/"
load(paste0(out_folder, "Writeup11_rosmap-initial.RData"))

seurat_obj
summary(seurat_obj@meta.data)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)

seurat_obj
summary(seurat_obj@meta.data)

# replace the misisng PMI with the median
vec <- seurat_obj$pmi
vec[is.na(vec)] <- stats::median(vec, na.rm = TRUE)
seurat_obj$pmi <- vec

zz <- SeuratObject::LayerData(seurat_obj,
                              assay = "RNA", 
                              layer = "counts")
SeuratObject::LayerData(seurat_obj,
                        assay = "RNA", 
                        layer = "data") <- zz

# remove "data" slot
tmp <- SeuratObject::LayerData(object = seurat_obj, 
                               assay = "RNA", 
                               layer = "counts") 
head(tmp@x)

# convert everything back into characters in preparation for the conversion
variable_names <- c("subject", "batch", "ADdiag3types", "seurat_clusters", "celltype", 
                    "Pt_ID", "Study", "race", "spanish", "apoe_genotype", 
                    "sex", "APOEe4_status", "ADpath", "brainRegion")
for(variable in variable_names){
  seurat_obj@meta.data[,variable] <- as.character(seurat_obj@meta.data[,variable])
}
summary(seurat_obj@meta.data)

# https://github.com/mojaveazure/seurat-disk/issues/27
# to appease the conversion, you need variable features
Seurat::VariableFeatures(seurat_obj) <- SeuratObject::Features(seurat_obj)

# https://github.com/satijalab/seurat/issues/8220#issuecomment-1871874649
seurat_obj[["RNA3"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
Seurat::DefaultAssay(seurat_obj) <- "RNA3"
seurat_obj[["RNA"]] <- NULL
seurat_obj <- SeuratObject::RenameAssays(object = seurat_obj, RNA3 = 'RNA')

SeuratDisk::SaveH5Seurat(seurat_obj, 
                         filename = paste0(out_folder, "Writeup11_rosmap-initial.h5Seurat"))
SeuratDisk::Convert(paste0(out_folder, "Writeup11_rosmap-initial.h5Seurat"),
                    dest = "h5ad", 
                    misc = FALSE)

print("Done! :)")