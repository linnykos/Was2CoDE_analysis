rm(list=ls())
library(Seurat)
library(eSVD2)
set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/naive-preprocess.RData")

#########

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
categorical_vars <- c("Sex", "SeqBatch", "Race", "CognitiveStatus")
numerical_vars <- c("PMI", "BrainPh", "FreshBrainWeight", "coded_Age")

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

na_vars <- c("PMI", "BrainPh", "Race")
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

###########

gene_vec <- Seurat::VariableFeatures(ss_data_norm[["RNA"]])
ss_data_norm <- subset(ss_data_norm, features = gene_vec)

mat <- Matrix::t(SeuratObject::LayerData(ss_data_norm, 
                                         assay = "RNA", 
                                         layer = "counts"))

covariate_dat <- ss_data_norm@meta.data[,c("Pt_ID", categorical_vars, numerical_vars)]
covariate_df <- data.frame(covariate_dat)
for(variable in setdiff(c("Pt_ID", categorical_vars), "CognitiveStatus")){
  covariate_df[,variable] <- factor(covariate_df[,variable], levels = names(sort(table(covariate_df[,variable]), decreasing = T)))
}
covariate_df[,"CognitiveStatus"] <- factor(covariate_df[,"CognitiveStatus"], levels = c("No_dementia", "Dementia"))
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = numerical_vars)

############

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariates[,-grep("Pt_ID", colnames(covariates))],
                                    case_control_variable = "CognitiveStatus_Dementia",
                                    bool_intercept = T,
                                    k = 30,
                                    lambda = 0.1,
                                    metadata_case_control = covariates[,"CognitiveStatus_Dementia"],
                                    metadata_individual = covariate_df[,"Pt_ID"],
                                    verbose = 1)
time_end1 <- Sys.time()

omitted_variables <- colnames(eSVD_obj$covariates)[grep("SeqBatch", colnames(eSVD_obj$covariates))]
eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = c("Log_UMI", omitted_variables)
)

print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = setdiff(colnames(eSVD_obj$covariates), "diagnosis_ASD"),
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_First",
                             fit_previous = "fit_Init")
time_end2 <- Sys.time()

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_First",
  omitted_variables = c("Log_UMI", omitted_variables)
)

print("Second fit")
time_start3 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = NULL,
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_Second",
                             fit_previous = "fit_First")
time_end3 <- Sys.time()

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Second",
  omitted_variables = omitted_variables
)

print("Nuisance estimation")
time_start4 <- Sys.time()
eSVD_obj <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                      bool_covariates_as_library = T,
                                      bool_library_includes_interept = T,
                                      bool_use_log = F,
                                      verbose = 1)
time_end4 <- Sys.time()

eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      alpha_max = 2*max(mat@x),
                                      bool_covariates_as_library = T,
                                      bool_stabilize_underdispersion = T,
                                      library_min = 0.1,
                                      pseudocount = 0)

time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           verbose = 1)
eSVD_obj <- eSVD2:::compute_pvalue(input_obj = eSVD_obj)
time_end5 <- Sys.time()


#########

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/preprocess/naive-preprocess.RData.",
              "Applying eSVD2.")

save(date_of_run, session_info, note,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../../out/tati/Writeup2/Writeup2_esvd.RData")

print("Done! :)")
