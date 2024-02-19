rm(list=ls())
library(openxlsx)

metadata <- openxlsx::read.xlsx("~/kzlinlab/data/sea-ad/sea-ad_cohort_donor_metadata_082222.xlsx")
rownames(metadata) <- metadata$Donor.ID
harmonized_scores <- openxlsx::read.xlsx("~/kzlinlab/data/sea-ad/joey-harmonized-cognitive/sea_ad_cognitive_slopes_79.xlsx")
rownames(harmonized_scores) <- harmonized_scores$donor_name
neuropath <- read.csv("~/kzlinlab/data/sea-ad/sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv")
rownames(neuropath) <- neuropath$Donor.ID
path_idx <- grep("^percent.*area_Grey", colnames(neuropath))
colnames(neuropath)[path_idx]
neuropath <- neuropath[,path_idx]

included_idx <- which(metadata$Donor.ID %in% harmonized_scores$donor_name)
metadata <- metadata[included_idx,]

ad_idx <- intersect(
  which(metadata$Braak %in% c("Braak V", "Braak VI")),
  which(metadata$CERAD.score %in% c("Moderate", "Frequent"))
)
metadata <- metadata[ad_idx,]
harmonized_scores <- harmonized_scores[rownames(metadata),]
neuropath <- neuropath[rownames(metadata),]

# create the cognitive covariate
cognition_vec <- rep("AD", nrow(harmonized_scores))
names(cognition_vec) <- metadata$Donor.ID
med_val <- stats::median(harmonized_scores$slope_zmem0)
high_donor <- rownames(harmonized_scores)[harmonized_scores$slope_zmem0 >= med_val]
cognition_vec[high_donor] <- "Resilient"

########################

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
set.seed(10)

# process the seurat_obj$donor_id
donor_vec <- seurat_obj$donor_id
donor_vec <- sapply(donor_vec, function(x){
  paste0(substr(x, start = 0, stop = 3), ".", 
         substr(x, start = 4, stop = 5), ".",
         substr(x, start = 6, stop = 8))
})
seurat_obj$donor_id <- donor_vec

resilient_vec <- rep(NA, length(SeuratObject::Cells(seurat_obj)))
resilient_vec[seurat_obj$donor_id %in% names(cognition_vec)[which(cognition_vec == "AD")]] <- "AD"
resilient_vec[seurat_obj$donor_id %in% names(cognition_vec)[which(cognition_vec == "Resilient")]] <- "Resilient"
seurat_obj$resilient <- resilient_vec

# first subset only the individuals in our analysis
keep_vec <- rep(FALSE, length(SeuratObject::Cells(seurat_obj)))
keep_vec[which(seurat_obj$donor_id %in% names(cognition_vec))] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# construct the age vector
age_vec <- sapply(seurat_obj$development_stage, function(x){
  substr(x, start = 0, stop = 2)
})
age_vec[which(seurat_obj$development_stage == "80yearoldandoverhumanstage")] <- "90"
age_vec <- as.numeric(age_vec)
seurat_obj$AgeAtDeath <- age_vec

seurat_obj$Lewy.body.disease.pathology <- seurat_obj@meta.data[,"Lewy body disease pathology"]
seurat_obj$APOE4.status <- seurat_obj@meta.data[,"APOE4 status"]

# adjust the covariates
# need to convert to numeric
categorical_vars <- c("sex", "self_reported_ethnicity", "YearsOfEducation", "APOE4.status", "Lewy.body.disease.pathology", "resilient")
numerical_vars <- c("FractionMitochrondrialUMIs", "AgeAtDeath")

for(variable in categorical_vars){
  seurat_obj@meta.data[,variable] <- factor(seurat_obj@meta.data[,variable])
}

# now start the actualy processing
Seurat::DefaultAssay(seurat_obj) <- "RNA"
mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                         assay = "RNA", 
                                         features = Seurat::VariableFeatures(seurat_obj),
                                         layer = "counts"))

covariate_dat <- seurat_obj@meta.data[,c("donor_id", categorical_vars, numerical_vars)]
covariate_df <- data.frame(covariate_dat)
for(variable in setdiff(c("donor_id", categorical_vars), "resilient")){
  covariate_df[,variable] <- factor(covariate_df[,variable], levels = names(sort(table(covariate_df[,variable]), decreasing = T)))
}
covariate_df[,"resilient"] <- factor(covariate_df[,"resilient"], levels = c("AD", "Resilient"))
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = numerical_vars)

############

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariates[,-grep("donor_id", colnames(covariates))],
                                    case_control_variable = "resilient_Resilient",
                                    bool_intercept = T,
                                    k = 30,
                                    lambda = 0.1,
                                    metadata_case_control = covariates[,"resilient_Resilient"],
                                    metadata_individual = covariate_df[,"donor_id"],
                                    verbose = 1)
time_end1 <- Sys.time()

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
                             offset_variables = setdiff(colnames(eSVD_obj$covariates), "resilient_Resilient"),
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
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData.",
              "Applying eSVD2. Using only subjects with AD pathology.")

save(date_of_run, session_info, note,
     eSVD_obj,
     harmonized_scores, neuropath,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")

print("Done! :)")
