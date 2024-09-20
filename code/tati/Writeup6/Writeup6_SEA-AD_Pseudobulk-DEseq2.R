rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData") 

categorical_vars <- c("ADNC", "sex", "assay", "self_reported_ethnicity", "APOE4status")
numerical_vars <- c("Ageatdeath", "PMI")

# Now extract the data
zz <- seurat_obj@meta.data[, c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)


###########
# gene_vec <- Seurat::VariableFeatures(seurat_obj[["RNA"]])
# seurat_obj <- subset(seurat_obj, features = gene_vec)

########

subj_vec <- seurat_obj$donor_id

Seurat::DefaultAssay(seurat_obj) <- "RNA"
mat <- SeuratObject::LayerData(seurat_obj, 
                               assay = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(seurat_obj))
mat <- as.matrix(mat)
# Prepare metadata
metadata <- seurat_obj@meta.data[,c(categorical_vars, numerical_vars)]
for(var in categorical_vars){
  metadata[,var] <- as.factor(metadata[,var])
}

# Create a unique identifier for each donor_id x assay pair
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
  mat_pseudobulk[,superstring] <- Matrix::rowSums(mat[,idx])
  
  for(vr in categorical_vars){
    metadata_pseudobulk[superstring, vr] <- as.character(unique(metadata[idx,vr]))
  }
  for(vr in numerical_vars){
    metadata_pseudobulk[superstring, vr] <- mean(metadata[idx,vr])
  }
}

metadata_pseudobulk[,"ADNC"] <- factor(metadata_pseudobulk[,"ADNC"])
metadata_pseudobulk[,"ADNC"] <- relevel(metadata_pseudobulk[,"ADNC"], ref = "Control")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ ADNC + sex + PMI + assay + Ageatdeath + APOE4status + self_reported_ethnicity)

dds <- DESeq2::DESeq(dds) #
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="ADNC_Case_vs_Control")
deseq2_res$padj_custom <- stats::p.adjust(deseq2_res[,"pvalue"], method = "BH")

save(deseq2_res,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_Pseudobulk-DEseq2.RData")
