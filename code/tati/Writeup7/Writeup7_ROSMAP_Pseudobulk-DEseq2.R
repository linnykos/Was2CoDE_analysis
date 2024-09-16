rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
library(DESeq2)

set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup11/Writeup11_rosmap_scVI-postprocessed.RData") 

# table(seurat_obj$Pt_ID, seurat_obj$APOEe4_status)
categorical_vars <- c("ADpath", "sex", "race", "batch", "APOEe4_status")
numerical_vars <- c("age_death", "pmi")

zz <- seurat_obj@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)

###########
# gene_vec <- Seurat::VariableFeatures(seurat_obj[["RNA"]])
# seurat_obj <- subset(seurat_obj, features = gene_vec)

########

subj_vec <- seurat_obj$Pt_ID

Seurat::Defaultbatch(seurat_obj) <- "RNA"
mat <- SeuratObject::LayerData(seurat_obj, 
                               batch = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(seurat_obj))
mat <- as.matrix(mat)
# Prepare metadata
metadata <- seurat_obj@meta.data[,c(categorical_vars, numerical_vars)]
for(var in categorical_vars){
  metadata[,var] <- as.factor(metadata[,var])
}

# Create a unique identifier for each Pt_ID x batch pair
superstring_vec <- sapply(1:nrow(metadata), function(i){
  paste0(subj_vec[i], "-", sapply(categorical_vars, function(var){
    as.character(metadata[i,var])
  }), collapse = "-")
})
unique_superstring <- unique(superstring_vec)
num_uniq <- length(unique_superstring)

# Create pseudobulk matrix and associated metadata
mat_pseudobulk <- matrix(NA, nrow = nrow(mat), ncol = num_uniq)
rownames(mat_pseudobulk) <- rownames(mat)
colnames(mat_pseudobulk) <- unique_superstring
metadata_pseudobulk <- as.data.frame(matrix(NA, nrow = num_uniq, ncol = ncol(metadata)))
colnames(metadata_pseudobulk) <- colnames(metadata)
rownames(metadata_pseudobulk) <- unique_superstring


for(superstring in unique_superstring){
  idx <- which(superstring_vec == superstring)
  
  if(length(idx) == 1){
    # If only one column, directly assign it
    mat_pseudobulk[,superstring] <- mat[, idx]
  } else {
    # Otherwise, calculate rowSums
    mat_pseudobulk[,superstring] <- Matrix::rowSums(mat[, idx])
  }
  
  for(vr in categorical_vars){
    metadata_pseudobulk[superstring, vr] <- as.character(unique(metadata[idx, vr]))
  }
  
  for(vr in numerical_vars){
    metadata_pseudobulk[superstring, vr] <- mean(metadata[idx, vr])
  }
}


metadata_pseudobulk[,"ADpath"] <- factor(metadata_pseudobulk[,"ADpath"])
metadata_pseudobulk[,"ADpath"] <- relevel(metadata_pseudobulk[,"ADpath"], ref = "no")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ ADpath + sex + pmi + batch + age_death + APOEe4_status + race)

dds <- DESeq2::DESeq(dds) #this is the function that takes a while
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)
#resultsNames(dds)
deseq2_res <- DESeq2::results(dds, name="ADpath_yes_vs_no")
deseq2_res$padj_custom <- stats::p.adjust(deseq2_res[,"pvalue"], method = "BH")

save(deseq2_res,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup7/Writeup7_ROSMAP_Pseudobulk-DEseq2.RData")
