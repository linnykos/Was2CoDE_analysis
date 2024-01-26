rm(list=ls())
library("Seurat")

seurat_obj <- readRDS("~/kzlinlab/data/sea-ad/microglia-pvm_dpc.rds")

#########

seurat_obj$CognitiveStatus <- seurat_obj$`Cognitive status`  
seurat_obj$`Cognitive status` <- NULL
seurat_obj$AgeAtDeath <- seurat_obj$`Age at death`  
seurat_obj$`Age at death` <- NULL
seurat_obj$YearsOfEducation <- seurat_obj$`Years of education`  
seurat_obj$`Years of education` <- NULL
seurat_obj$FractionMitochrondrialUMIs <- seurat_obj$`Fraction mitochrondrial UMIs`  
seurat_obj$`Fraction mitochrondrial UMIs` <- NULL

tmp <- seurat_obj$CognitiveStatus
seurat_obj$CognitiveStatus <- factor(gsub(
  pattern = " ",
  replacement = "_",
  x = tmp
))

#########

# adjust the covariates
# need to convert to numeric
categorical_vars <- c("sex", "self_reported_ethnicity", "AgeAtDeath", "YearsOfEducation", "CognitiveStatus")
numerical_vars <- c("FractionMitochrondrialUMIs")

for(variable in categorical_vars){
  seurat_obj@meta.data[,variable] <- factor(seurat_obj@meta.data[,variable])
}

###########

gene_vec <- Seurat::VariableFeatures(seurat_obj[["RNA"]])
seurat_obj <- subset(seurat_obj, features = gene_vec)

# remove the reference cells
keep_vec <- rep(FALSE, ncol(seurat_obj))
keep_vec[seurat_obj$CognitiveStatus != "Reference"] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                         assay = "RNA", 
                                         layer = "counts"))

covariate_dat <- seurat_obj@meta.data[,c("donor_id", categorical_vars, numerical_vars)]
covariate_df <- data.frame(covariate_dat)
for(variable in setdiff(c("donor_id", categorical_vars), "CognitiveStatus")){
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
                                    covariates = covariates[,-grep("donor_id", colnames(covariates))],
                                    case_control_variable = "CognitiveStatus_Dementia",
                                    bool_intercept = T,
                                    k = 30,
                                    lambda = 0.1,
                                    metadata_case_control = covariates[,"CognitiveStatus_Dementia"],
                                    metadata_individual = covariate_df[,"donor_id"],
                                    verbose = 1)
time_end1 <- Sys.time()


save(eSVD_obj,
     file = "../../../../../out/kevin/Writeup1/Writeup1_esvd.RData")

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = "Log_UMI"
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
  omitted_variables = "Log_UMI"
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
  omitted_variables = NULL
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
note <- paste("Working from ~/lab/projects/subject-de/out/kevin/Writeup1/naive-preprocess.RData.",
              "Applying eSVD2.")

save(date_of_run, session_info, note,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../../out/kevin/Writeup1/Writeup1_esvd.RData")

print("Done! :)")

