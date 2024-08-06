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

categorical_vars <- c("Sex", "SeqBatch", "Race", "Study_Designation")
numerical_vars <- c("PMI", "BrainPh", "FreshBrainWeight", "coded_Age","genotype_APOE")

tmp <- as.character(ss_data_norm$coded_Age)
tmp[which(tmp == "90+")] <- "90"
ss_data_norm$coded_Age <- tmp

for(variable in categorical_vars){
  ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
}
for(variable in numerical_vars){
  ss_data_norm@meta.data[,variable] <- as.numeric(as.character(ss_data_norm@meta.data[,variable]))
}

# need to fill in NAs in PMI, BrainPh, Race
tab_list <- lapply(categorical_vars, function(variable){
  table(ss_data_norm$Pt_ID, ss_data_norm@meta.data[,variable])
})
names(tab_list) <- categorical_vars

na_vars <- c("PMI", "BrainPh", "Race", "genotype_APOE")
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

##########################

neb_data <- nebula::scToNeb(obj = ss_data_norm,
                            assay = "RNA",
                            id = "Pt_ID",
                            pred = c("Study_Designation", "Sex", "PMI", "SeqBatch", "coded_Age", "genotype_APOE"),
                            offset = "nCount_RNA")
df <- model.matrix( ~ Study_Designation + Sex + PMI + SeqBatch + coded_Age + genotype_APOE,
                    data = neb_data$pred)
start_time <- Sys.time()
nebula_res <- nebula::nebula(count = neb_data$count,
                             id = neb_data$id,
                             pred = df,
                             offset = neb_data$offset,
                             model = "NBGMM",
                             verbose = TRUE)
end_time <- Sys.time()

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/preprocess/naive-preprocess.RData.",
              "Applying NEBULA.")

save(nebula_res, 
     date_of_run, session_info, note,
     start_time, end_time,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup3/Writeup3_NEBULA.RData")

print("Done! :)")