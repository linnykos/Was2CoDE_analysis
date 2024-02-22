rm(list=ls())
library(Seurat)
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
colname_vec <- colnames(neuropath)
colname_vec <- sapply(colname_vec, function(val){
  tmp <- strsplit(val, split = "\\.")[[1]][1:2]
  paste0(tmp, collapse = ".")
})
colnames(neuropath) <- colname_vec

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
seurat_obj@meta.data[,"Lewy body disease pathology"] <- NULL
seurat_obj@meta.data[,"APOE4 status"] <- NULL

#############################

# do some cleanup for scVI later
categorical_vars <- c("donor_id", "sex", "self_reported_ethnicity", "YearsOfEducation", "APOE4.status", "Lewy.body.disease.pathology", "resilient")
numerical_vars <- c("FractionMitochrondrialUMIs", "AgeAtDeath")
all_vars <- c(categorical_vars, numerical_vars)
metadata_vars <- colnames(seurat_obj@meta.data)
metadata_vars <- metadata_vars[!metadata_vars %in% all_vars]
for(variable in metadata_vars){
  seurat_obj@meta.data[,variable] <- NULL
}
for(variable in all_vars){
  vec <- seurat_obj@meta.data[,variable]
  if(is.character(vec)) vec <- factor(vec)
  if(is.factor(vec)) vec <- droplevels(vec)
  seurat_obj@meta.data[,variable] <- vec
}
colname_vec <- colnames(neuropath)
donor_vec <- seurat_obj$donor_id
for(variable in colname_vec){
  seurat_obj@meta.data[,variable] <- neuropath[donor_vec,variable]
}

seurat_diet <- Seurat::DietSeurat(
  seurat_obj,
  layers = "counts",
  features = NULL, 
  assays = "RNA",
  graphs = NULL, 
  misc = FALSE
)

options(Seurat.object.assay.version = "v3")
seurat_diet[["RNA"]] <- as(object = seurat_diet[["RNA"]], Class = "Assay")

#######

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Input from ~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData.",
              "Added neuropath data, did some reformatting, and changed to Seuratv3 since I'm going to",
              "use a older version of Seurat to convert this into an anndata object")

save(date_of_run, session_info, 
     seurat_diet,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess_Seuratv3.RData")
