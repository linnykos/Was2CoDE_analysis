rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/naive-preprocess.RData")

gene_vec <- Seurat::VariableFeatures(ss_data_norm[["RNA"]])
ss_data_norm <- subset(ss_data_norm, features = gene_vec)

# adjust the covariates
# need to convert to numeric
tmp <- paste0("ID_", as.character(ss_data_norm$Pt_ID))
ss_data_norm$Pt_ID <- factor(tmp)

categorical_vars <- c("Sex", "SeqBatch", "Race", "Study_Designation","genotype_APOE")
numerical_vars <- c("PMI","coded_Age")

tmp <- as.character(ss_data_norm$coded_Age)
tmp[which(tmp == "90+")] <- "90"
ss_data_norm$coded_Age <- tmp

for(variable in categorical_vars){
  ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
}
for(variable in numerical_vars){
  ss_data_norm@meta.data[,variable] <- as.numeric(as.character(ss_data_norm@meta.data[,variable]))
}

# need to fill in NAs in PMI,Race
tab_list <- lapply(categorical_vars, function(variable){
  table(ss_data_norm$Pt_ID, ss_data_norm@meta.data[,variable])
})
names(tab_list) <- categorical_vars

pmi_missing_values <- c(ID_1 = 5.08,
                        ID_10 = 3.33,
                        ID_3 = 5.08,
                        ID_15 = 6.97)
for(i in 1:length(pmi_missing_values)){
  subj_id <- names(pmi_missing_values)[i]
  idx <- which(ss_data_norm$Pt_ID == subj_id)
  ss_data_norm$PMI[idx] <- pmi_missing_values[i]
}

race_missing_values <- c(ID_5 = "White")
for(i in 1:length(race_missing_values)){
  subj_id <- names(race_missing_values)[i]
  idx <- which(ss_data_norm$Pt_ID == subj_id)
  ss_data_norm$Race[idx] <- race_missing_values[i]
}

zz <- ss_data_norm@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)



########

subj_vec <- ss_data_norm$Pt_ID

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
mat <- SeuratObject::LayerData(ss_data_norm, 
                               assay = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(ss_data_norm))
mat <- as.matrix(mat)
categorical_var <- c("Study_Designation", "Race","Pt_ID", "Sex", "SeqBatch", "genotype_APOE")
numerical_var <- c("PMI", "coded_Age")
metadata <- ss_data_norm@meta.data[,c(categorical_var, numerical_var)]
for(var in categorical_var){
  metadata[,var] <- as.factor(metadata[,var])
}

uniq_subj_vec <- unique(subj_vec)
num_subj <- length(uniq_subj_vec)
mat_pseudobulk <- matrix(NA, nrow = nrow(mat), ncol = num_subj)
rownames(mat_pseudobulk) <- rownames(mat)
colnames(mat_pseudobulk) <- uniq_subj_vec
metadata_pseudobulk <- as.data.frame(matrix(NA, nrow = num_subj, ncol = ncol(metadata)))
colnames(metadata_pseudobulk) <- colnames(metadata)
rownames(metadata_pseudobulk) <- uniq_subj_vec

for(subj in uniq_subj_vec){
  idx <- which(subj_vec == subj)
  mat_pseudobulk[,subj] <- Matrix::rowSums(mat[,idx])
  
  for(vr in categorical_var){
    metadata_pseudobulk[subj, vr] <- as.character(unique(metadata[idx,vr]))
  }
  for(vr in numerical_var){
    metadata_pseudobulk[subj, vr] <- mean(metadata[idx,vr])
  }
}

for(vr in categorical_var){
  metadata_pseudobulk[,vr] <- factor(metadata_pseudobulk[,vr])
}
for(vr in numerical_var){
  metadata_pseudobulk[,vr] <- scale(metadata_pseudobulk[,vr])
}

metadata_pseudobulk[,"Study_Designation"] <- relevel(metadata_pseudobulk[,"Study_Designation"], ref = "No_dementia")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ Study_Designation + Sex + PMI + SeqBatch + coded_Age + genotype_APOE + Race)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="Study_Designation_Dementia_vs_No_dementia")
deseq2_res$padj_custom <- stats::p.adjust(deseq2_res[,"pvalue"], method = "BH")

save(deseq2_res,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_Katie_Pseudobulk-DEseq2.RData")