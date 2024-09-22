rm(list=ls())
library(Seurat)

out_folder <- "~/Dropbox/Collaboration-and-People/tati/out/kevin/Writeup13/"
load(paste0(out_folder, "Writeup13_simulation_data.RData"))

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- Matrix::t(SeuratObject::LayerData(seurat_obj,
                                         layer = "counts",
                                         assay = "RNA"))
covariate_dat <- seurat_obj@meta.data[,c("cc", "age", "gender", "tobacco", "individual")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"cc"] <- as.factor(covariate_df[,"cc"])
covariate_df[,"tobacco"] <- factor(covariate_df[,"tobacco"], levels = names(sort(table(covariate_df[,"tobacco"]), decreasing = T)))
covariate_df[,"gender"] <- factor(covariate_df[,"gender"], levels = names(sort(table(covariate_df[,"gender"]), decreasing = T)))
covariate_df[,"individual"] <- factor(covariate_df[,"individual"], levels = names(sort(table(covariate_df[,"individual"]), decreasing = T)))
covariates <- eSVD2::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("age"))

print("Initialization")
eSVD_obj <- eSVD2::initialize_esvd(dat = mat,
                                    covariates = covariates[,-grep("individual", colnames(covariates))],
                                    case_control_variable = "cc_1",
                                    bool_intercept = T,
                                    k = 10,
                                    lambda = 0.1,
                                    metadata_case_control = covariates[,"cc_1"],
                                    metadata_individual = covariate_df[,"individual"],
                                    verbose = 1)

eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = "Log_UMI"
)

print("First fit")
eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = setdiff(colnames(eSVD_obj$covariates), "cc_1"),
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_First",
                             fit_previous = "fit_Init")

eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_First",
  omitted_variables = "Log_UMI"
)

print("Second fit")
eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = NULL,
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_Second",
                             fit_previous = "fit_First")

eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Second",
  omitted_variables = numeric(0)
)

print("Nuisance estimation")
eSVD_obj <- eSVD2::estimate_nuisance(input_obj = eSVD_obj,
                                      bool_covariates_as_library = TRUE,
                                      bool_library_includes_interept = TRUE,
                                      bool_use_log = FALSE,
                                      verbose = 1)

eSVD_obj <- eSVD2::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = FALSE,
                                      alpha_max = NULL,
                                      bool_covariates_as_library = TRUE,
                                      bool_stabilize_underdispersion = TRUE,
                                      library_min = 1,
                                      pseudocount = 1)
eSVD_obj <- eSVD2::compute_test_statistic(input_obj = eSVD_obj,
                                           verbose = 1)
eSVD_obj <- eSVD2::compute_pvalue(input_obj = eSVD_obj)

##############

save(eSVD_obj,
     seurat_obj,
     case_individuals,
     control_individuals,
     covariates,
     gene_labeling,
     gene_labeling2,
     gene_library_vec,
     individual_vec,
     nuisance_vec,
     nat_mat,
     obs_mat,
     true_fdr_vec,
     true_logpvalue_vec,
     true_null_mean,
     true_null_sd,
     true_teststat_vec,
     x_mat,
     y_mat,
     z_mat,
     date_of_run, session_info,
     file = paste0(out_folder, "Writeup13_esvd.RData"))

