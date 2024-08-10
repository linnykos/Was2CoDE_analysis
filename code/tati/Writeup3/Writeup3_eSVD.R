rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
library(eSVD2)
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
###########

gene_vec <- Seurat::VariableFeatures(ss_data_norm[["RNA"]])
ss_data_norm <- subset(ss_data_norm, features = gene_vec)

mat <- Matrix::t(SeuratObject::LayerData(ss_data_norm, 
                                         assay = "RNA", 
                                         layer = "counts"))

covariate_dat <- ss_data_norm@meta.data[,c("Pt_ID", categorical_vars, numerical_vars)]
covariate_df <- data.frame(covariate_dat)
for(variable in setdiff(c("Pt_ID", categorical_vars), "Study_Designation")){
  covariate_df[,variable] <- factor(covariate_df[,variable], levels = names(sort(table(covariate_df[,variable]), decreasing = T)))
}
covariate_df[,"Study_Designation"] <- factor(covariate_df[,"Study_Designation"], levels = c("Ctrl", "AD"))
covariates <- eSVD2::format_covariates(dat = mat,
                                       covariate_df = covariate_df,
                                       rescale_numeric_variables = numerical_vars)

############

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2::initialize_esvd(dat = mat,
                                   covariates = covariates[,-grep("Pt_ID", colnames(covariates))],
                                   case_control_variable = "Study_Designation_AD",
                                   bool_intercept = T,
                                   k = 30,
                                   lambda = 0.1,
                                   metadata_case_control = covariates[,"Study_Designation_AD"],
                                   metadata_individual = covariate_df[,"Pt_ID"],
                                   verbose = 1)
time_end1 <- Sys.time()

omitted_variables <- colnames(eSVD_obj$covariates)[grep("SeqBatch", colnames(eSVD_obj$covariates))]
eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = c("Log_UMI", omitted_variables)
)

print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                            l2pen = 0.1,
                            max_iter = 100,
                            offset_variables = setdiff(colnames(eSVD_obj$covariates), "diagnosis_ASD"),
                            tol = 1e-6,
                            verbose = 1,
                            fit_name = "fit_First",
                            fit_previous = "fit_Init")
time_end2 <- Sys.time()

eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_First",
  omitted_variables = c("Log_UMI", omitted_variables)
)

print("Second fit")
time_start3 <- Sys.time()
eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                            l2pen = 0.1,
                            max_iter = 100,
                            offset_variables = NULL,
                            tol = 1e-6,
                            verbose = 1,
                            fit_name = "fit_Second",
                            fit_previous = "fit_First")
time_end3 <- Sys.time()

eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Second",
  omitted_variables = omitted_variables
)

print("Nuisance estimation")
time_start4 <- Sys.time()
eSVD_obj <- eSVD2::estimate_nuisance(input_obj = eSVD_obj,
                                     bool_covariates_as_library = T,
                                     bool_library_includes_interept = T,
                                     bool_use_log = F,
                                     verbose = 1)
time_end4 <- Sys.time()

eSVD_obj <- eSVD2::compute_posterior(input_obj = eSVD_obj,
                                     bool_adjust_covariates = F,
                                     alpha_max = 2*max(mat@x),
                                     bool_covariates_as_library = T,
                                     bool_stabilize_underdispersion = T,
                                     library_min = 0.1,
                                     pseudocount = 0)

time_start5 <- Sys.time()
eSVD_obj <- eSVD2::compute_test_statistic(input_obj = eSVD_obj,
                                          verbose = 1)
eSVD_obj <- eSVD2::compute_pvalue(input_obj = eSVD_obj)
time_end5 <- Sys.time()


#########

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData.",
              "Applying eSVD2.")

save(date_of_run, session_info, note,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup3/Writeup3_prater_esvd.RData")

print("Done! :)")
