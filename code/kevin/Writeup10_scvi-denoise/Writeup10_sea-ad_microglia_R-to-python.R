rm(list=ls())
library("Seurat")

seurat_obj <- readRDS("~/kzlinlab/data/sea-ad/microglia-pvm_dpc.rds")
seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)

#########

# remove unnecessary objects
Seurat::DefaultAssay(seurat_obj) <- "RNA"
SeuratObject::LayerData(object = seurat_obj, 
                        assay = "RNA", 
                        layer = "data") <- NULL

# go through all the column names and remove all punctuations (aside from "_" and ".")
colnames_vec <- colnames(seurat_obj@meta.data)
colnames_vec <- gsub(pattern = "[^[:alnum:]._ ]", replacement = "", x = colnames_vec)
colnames_vec <- gsub(pattern = " ", replacement = "", x = colnames_vec)
colnames(seurat_obj@meta.data) <- colnames_vec

# go through each column. If it's a numeric, ignore. If it's a factor or character, remove all punctuations (aside from "_" and ".")
for(j in 1:ncol(seurat_obj@meta.data)){
  tmp <- seurat_obj@meta.data[,j]
  if(is.numeric(tmp)) next()
  if(is.factor(tmp)) tmp <- as.character(tmp)
  tmp <- gsub(pattern = "[^[:alnum:]._ ]", replacement = "", x = tmp)
  tmp <- gsub(pattern = " ", replacement = "", x = tmp)
  seurat_obj@meta.data[,j] <- tmp
}

# Fix the age vector
age_vec <- sapply(seurat_obj$development_stage, function(x){
  substr(x, start = 0, stop = 2)
})
age_vec[which(seurat_obj$development_stage == "80yearoldandoverhumanstage")] <- "90"
age_vec <- as.numeric(age_vec)
seurat_obj$Ageatdeath <- age_vec

# remove the donors that are "Reference"
keep_vec <- rep(TRUE, length(Seurat::Cells(seurat_obj)))
keep_vec[which(seurat_obj$ADNC == "Reference")] <- FALSE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# rename the current embeddings
seurat_obj[["original_scVI"]] <- seurat_obj[["scVI"]]
seurat_obj[["original_umap"]] <- seurat_obj[["umap"]]
seurat_obj[["scVI"]] <- NULL
seurat_obj[["umap"]] <- NULL

# fix PMI
vec <- seurat_obj$PMI
vec <- sapply(strsplit(vec, split = "to"), function(x){x[1]})
seurat_obj$PMI <- as.numeric(vec)

# checking
categorical_vars <- c("sex", "assay", "self_reported_ethnicity", "ADNC", "APOE4status", "donor_id")
numerical_vars <- c("PMI","Ageatdeath")

zz <- seurat_obj@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(sapply(zz[,numerical_vars], is.numeric))
stopifnot(sapply(zz[,categorical_vars], is.character))
stopifnot(!any(is.na(zz)))
stopifnot(!any(zz == "NA"))
stopifnot(all(!sapply(1:ncol(zz), is.factor)))
summary(zz)
head(zz)

###################### # Saving as R

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Working from ~/kzlinlab/data/sea-ad/microglia-pvm_dpc.rds.")

save(date_of_run, session_info, note,
     seurat_obj,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_cleaned.RData")

###################### # Saving as python

SeuratDisk::SaveH5Seurat(seurat_obj, 
                         filename = "~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_cleaned.h5Seurat")
SeuratDisk::Convert("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_cleaned.h5Seurat", 
                    dest = "h5ad", 
                    misc = FALSE)
