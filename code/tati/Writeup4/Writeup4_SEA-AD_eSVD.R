rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
library(eSVD2)
set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")

# construct the age vector
age_vec <- sapply(seurat_obj$development_stage, function(x){
  substr(x, start = 0, stop = 2)
})

age_vec[which(seurat_obj$development_stage == "80yearoldandoverhumanstage")] <- "90"
age_vec <- as.numeric(age_vec)
seurat_obj$AgeAtDeath <- age_vec


seurat_obj$Lewy.body.disease.pathology <- seurat_obj@meta.data[,"Lewy body disease pathology"]
seurat_obj$APOE4.status <- seurat_obj@meta.data[,"APOE4 status"]
#colnames(seurat_obj@meta.data)

# gene_vec <- Seurat::VariableFeatures(seurat_obj@meta.data[["RNA"]])
# seurat_obj@meta.data <- subset(seurat_obj@meta.data, features = gene_vec)

tmp <- paste0("ID_", as.character(seurat_obj@meta.data$donor_id))
seurat_obj$donor_id <- factor(tmp)
# Remove "Reference" donors using subset function
seurat_obj <- subset(seurat_obj, subset = ADNC != "Reference")

# Reclassify ADNC levels
seurat_obj$ADNC <- with(seurat_obj@meta.data, 
                        ifelse(ADNC %in% c("NotAD", "Low"), "Control", 
                               ifelse(ADNC %in% c("Intermediate", "High"), "Case", ADNC)))

# table(seurat_obj$donor_id, seurat_obj$ADNC)

# Convert APOE4 status to binary (0 for "N", 1 for "Y")
seurat_obj$APOE4_status <- ifelse(seurat_obj@meta.data$APOE4.status == "Y", 1, 
                                  ifelse(seurat_obj@meta.data$APOE4.status == "N", 0, NA))

# table(seurat_obj$donor_id, seurat_obj$APOE4_status)
# Convert PMI Categories to (Numeric or) Factor
seurat_obj$PMI <- factor(seurat_obj$PMI, levels = c("32to59hours", "59to87hours", "87to114hours"))

categorical_vars <- c("ADNC", "sex", "assay", "self_reported_ethnicity","APOE4_status","PMI")
numerical_vars <- c("AgeAtDeath")

tmp <- as.character(seurat_obj$AgeAtDeath)
tmp[which(tmp == "90+")] <- "90"
seurat_obj@meta.data$AgeAtDeath <- tmp


for (variable in categorical_vars) {
  seurat_obj@meta.data[, variable] <- factor(seurat_obj@meta.data[, variable])
}

for (variable in numerical_vars) {
  seurat_obj@meta.data[, variable] <- as.numeric(as.character(seurat_obj@meta.data[, variable]))
}

zz <- seurat_obj@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)

###########
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
gene_vec <- Seurat::VariableFeatures(seurat_obj[["RNA"]])
seurat_obj <- subset(seurat_obj, features = gene_vec)

mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                         assay = "RNA", 
                                         layer = "counts"))

covariate_dat <- seurat_obj@meta.data[,c("donor_id", categorical_vars, numerical_vars)]
covariate_df <- data.frame(covariate_dat)

# table(seurat_obj$ADNC)
covariate_df$ADNC <- factor(seurat_obj$ADNC, levels = c("Control", "Case"))
for(variable in setdiff(c("donor_id", categorical_vars), "ADNC")){
  covariate_df[,variable] <- factor(covariate_df[,variable], levels = names(sort(table(covariate_df[,variable]), decreasing = TRUE)))
}
# colSums(is.na(covariate_df))
covariates <- eSVD2::format_covariates(dat = mat,
                                       covariate_df = covariate_df,
                                       rescale_numeric_variables = numerical_vars)

# colnames(covariates) <- gsub(",", "_", colnames(covariates))

# print(colnames(covariates))
############

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2::initialize_esvd(dat = mat,
                                   covariates = covariates[,-grep("donor_id", colnames(covariates))],
                                   case_control_variable = "ADNC_Case",
                                   bool_intercept = T,
                                   k = 30,
                                   lambda = 0.1,
                                   metadata_case_control = covariates[,"ADNC_Case"],
                                   metadata_individual = covariate_df[,"donor_id"],
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
                            offset_variables = setdiff(colnames(eSVD_obj$covariates), "ADNC_Case"),
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
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData,
              Applying eSVD2.")

save(date_of_run, session_info, note,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup4/Writeup4_SEA-AD_eSVD.RData")

print("Done! :)")
