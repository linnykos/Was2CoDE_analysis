rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

gene_vec <- Seurat::VariableFeatures(ss_data_norm[["RNA"]])
ss_data_norm <- subset(ss_data_norm, features = gene_vec)


# adjust the covariates
# need to convert to numeric
tmp <- paste0("ID_", as.character(ss_data_norm$Pt_ID))
ss_data_norm$Pt_ID <- factor(tmp)

categorical_vars <- c("Sex", "SeqBatch", "Race", "Study_Designation","APOEe4_status")
numerical_vars <- c("PMI","coded_Age")

tmp <- as.character(ss_data_norm$coded_Age)
ss_data_norm$coded_Age <- tmp

# Convert categorical variables to factors and numerical variables to numeric
for(variable in categorical_vars){
  ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
}
for(variable in numerical_vars){
  ss_data_norm@meta.data[,variable] <- as.numeric(as.character(ss_data_norm@meta.data[,variable]))
}

# Check for NA s
zz <- ss_data_norm@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)


########

subj_vec <- ss_data_norm$Pt_ID
Seurat::DefaultAssay(ss_data_norm) <- "RNA"
# Extract the count matrix
mat <- SeuratObject::LayerData(ss_data_norm, 
                               assay = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(ss_data_norm))
mat <- as.matrix(mat)
# Prepare metadata
metadata <- ss_data_norm@meta.data[,c(categorical_vars, numerical_vars)]
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


metadata_pseudobulk[,"Study_Designation"] <- factor(metadata_pseudobulk[,"Study_Designation"]) # Convert 'Study_Designation' to a factor before releveling
metadata_pseudobulk[,"Study_Designation"] <- relevel(metadata_pseudobulk[,"Study_Designation"], ref = "Ctrl")

# Create the DESeqDataSet object and Run the DESeq2 differential expression analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ Study_Designation + Sex + PMI + SeqBatch + coded_Age + APOEe4_status + Race)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds) # Get the names of the coefficients
# Extract the p-values from the results and calculate quantiles
deseq2_pval <- DESeq2::results(dds)$pvalue
pval_quantiles <- stats::quantile(deseq2_pval, na.rm = T)
print(pval_quantiles)
# Extract the results for the specific comparison: AD vs Ctrl
deseq2_res <- DESeq2::results(dds, name="Study_Designation_AD_vs_Ctrl")
deseq2_res$padj_custom <- stats::p.adjust(deseq2_res[,"pvalue"], method = "BH")
head(deseq2_res)

save(deseq2_res,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_Pseudobulk-DEseq2.RData")
