rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
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

########

subj_vec <- seurat_obj$donor_id

Seurat::DefaultAssay(ss_data_norm) <- "RNA"
mat <- SeuratObject::LayerData(ss_data_norm, 
                               assay = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(seurat_obj))
mat <- as.matrix(mat)
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

metadata_pseudobulk[,"ADNC"] <- relevel(metadata_pseudobulk[,"ADNC"], ref = "Control")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ ADNC + sex + PMI + assay + AgeAtDeath + APOE4_status + self_reported_ethnicity)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="ADNC_Case_vs_Control")
deseq2_res$padj_custom <- stats::p.adjust(deseq2_res[,"pvalue"], method = "BH")

save(deseq2_res,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup4/Writeup4_SEA-AD_Pseudobulk-DEseq2.RData")
