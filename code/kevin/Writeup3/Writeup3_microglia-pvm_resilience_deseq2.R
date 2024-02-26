rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

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

#########

subj_vec <- seurat_obj$donor_id

mat <- SeuratObject::LayerData(seurat_obj, 
                               assay = "RNA", 
                               features = Seurat::VariableFeatures(seurat_obj),
                               layer = "counts")
categorical_var <- c("sex", "self_reported_ethnicity", "YearsOfEducation", "APOE4.status", "Lewy.body.disease.pathology", "resilient")
numerical_var <- c("FractionMitochrondrialUMIs", "AgeAtDeath")

metadata <- seurat_obj@meta.data[,c(categorical_var, numerical_var)]
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

metadata_pseudobulk[,"resilient"] <- relevel(metadata_pseudobulk[,"resilient"], ref = "Resilient")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ sex + self_reported_ethnicity + YearsOfEducation + APOE4.status + Lewy.body.disease.pathology + FractionMitochrondrialUMIs + AgeAtDeath + resilient)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_res <- DESeq2::results(dds, name = "resilient_AD_vs_Resilient")

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(deseq2_res,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_deseq2.RData")
