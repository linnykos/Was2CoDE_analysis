rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

########

tmp <- paste0("ID_", as.character(ss_data_norm$Pt_ID))
ss_data_norm$Pt_ID <- factor(tmp)
tmp <- ss_data_norm$CognitiveStatus  
ss_data_norm$CognitiveStatus <- factor(gsub(
  pattern = " ",
  replacement = "_",
  x = tmp
))

#########

# adjust the covariates
# need to convert to numeric
categorical_vars <- c("Sex", "SeqBatch", "CognitiveStatus")
numerical_vars <- c("PMI", "coded_Age")

tmp <- as.character(ss_data_norm$coded_Age)
tmp[which(tmp == "90+")] <- "90"
ss_data_norm$coded_Age <- tmp

for(variable in categorical_vars){
  ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
}
for(variable in numerical_vars){
  ss_data_norm@meta.data[,variable] <- as.numeric(as.character(ss_data_norm@meta.data[,variable]))
}

# need to fill in NAs in PMI
tab_list <- lapply(categorical_vars, function(variable){
  table(ss_data_norm$Pt_ID, ss_data_norm@meta.data[,variable])
})
names(tab_list) <- categorical_vars

na_vars <- c("PMI")
num_nas <- sapply(na_vars, function(variable){
  tab_mat <- table(ss_data_norm$Pt_ID, ss_data_norm@meta.data[,variable])
  length(which(rowSums(tab_mat) == 0))
})
na_vars <- na_vars[order(num_nas, decreasing = F)]

for(variable in na_vars){
  tab_list <- lapply(categorical_vars, function(var_tmp){
    table(ss_data_norm$Pt_ID, ss_data_norm@meta.data[,var_tmp])
  })
  names(tab_list) <- categorical_vars
  
  tab_mat <- table(ss_data_norm$Pt_ID, ss_data_norm@meta.data[,variable])
  subj_to_match <- rownames(tab_mat)[which(rowSums(tab_mat) == 0)]
  
  # find all the categorical variables that match exactly
  subj_paired <- lapply(subj_to_match, function(subj){
    categorical_var_to_match <- setdiff(categorical_vars, variable)
    possible_idx <- which(table(unlist(lapply(categorical_var_to_match, function(variable_tmp){
      tab_mat <- tab_list[[variable_tmp]]
      col_idx <- which(tab_mat[subj,] != 0)
      idx <- which(tab_mat[,col_idx] != 0)
      rownames(tab_mat)[idx]
    }))) == length(categorical_var_to_match))
    
    possible_subj <- rownames(tab_list[[1]])[possible_idx]
    possible_subj <- setdiff(possible_subj, subj)
    stopifnot(length(possible_subj) > 0)
    possible_subj
  })
  names(subj_paired) <- subj_to_match
  
  # among these subjects, find the ones that are closest in numerical variables
  num_subjs <- length(subj_to_match)
  subj_paired2 <- sapply(1:num_subjs, function(ii){
    subj <- subj_to_match[ii]
    possible_subj <- subj_paired[[ii]]
    
    numerical_var_to_match <- setdiff(numerical_vars, variable)
    numerical_mat <- sapply(numerical_vars, function(var_tmp){
      sapply(c(subj, possible_subj), function(subj_tmp){
        tmp_idx <- which(ss_data_norm$Pt_ID == subj_tmp)
        tmp <- unique(ss_data_norm@meta.data[tmp_idx,var_tmp])
        stopifnot(length(tmp) == 1)
        tmp
      })
    })
    numerical_mat <- numerical_mat[,apply(numerical_mat, 2, function(x){!any(is.na(x))}), drop = F]
    numerical_mat <- scale(numerical_mat)
    dist_mat <- as.matrix(stats::dist(numerical_mat))
    diag(dist_mat) <- Inf
    rownames(dist_mat)[which.min(dist_mat[1,])]
  })
  names(subj_paired2) <- subj_to_match
  
  # now fill in values
  subject_vec <- ss_data_norm$Pt_ID
  zz <- ss_data_norm@meta.data[,variable]
  for(ii in 1:num_subjs){
    dest_idx <- which(subject_vec == names(subj_paired2)[ii])
    source_idx <- which(subject_vec == subj_paired2[ii])
    zz[dest_idx] <- unique(zz[source_idx])[1]
  }
  stopifnot(!any(is.na(zz)))
  if(variable %in% categorical_vars) zz <- factor(zz)
  ss_data_norm@meta.data[,variable] <- zz
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
categorical_var <- c("CognitiveStatus", "Pt_ID", "Sex", "SeqBatch")
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

metadata_pseudobulk[,"CognitiveStatus"] <- relevel(metadata_pseudobulk[,"CognitiveStatus"], ref = "No_dementia")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ Sex + SeqBatch + PMI + coded_Age + CognitiveStatus)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="CognitiveStatus_Dementia_vs_No_dementia")
deseq2_res$padj_custom <- stats::p.adjust(deseq2_res[,"pvalue"], method = "BH")

save(deseq2_res,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup2/Writeup2_Katie_Pseudobulk-DEseq2.RData")