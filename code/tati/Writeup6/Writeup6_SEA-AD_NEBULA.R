rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
set.seed(10)

print("Starting!")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData") 

# construct the age vector
age_vec <- sapply(seurat_obj$development_stage, function(x){
  substr(x, start = 0, stop = 2)
})

age_vec[which(seurat_obj$development_stage == "80yearoldandoverhumanstage")] <- "90"
age_vec <- as.numeric(age_vec)
seurat_obj$Ageatdeath <- age_vec


seurat_obj$Lewy.body.disease.pathology <- seurat_obj@meta.data[,"Lewybodydiseasepathology"]
seurat_obj$APOE4.status <- seurat_obj@meta.data[,"APOE4status"]
#colnames(seurat_obj@meta.data)

# gene_vec <- Seurat::VariableFeatures(seurat_obj@meta.data[["RNA"]])

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

categorical_vars <- c("ADNC", "sex", "assay", "self_reported_ethnicity","APOE4_status")
numerical_vars <- c("Ageatdeath","PMI")

tmp <- as.character(seurat_obj$Ageatdeath)
seurat_obj@meta.data$Ageatdeath <- tmp


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
gene_vec <- Seurat::VariableFeatures(seurat_obj[["RNA"]])
seurat_obj <- subset(seurat_obj, features = gene_vec)

##########################

neb_data <- nebula::scToNeb(obj = seurat_obj,
                            assay = "RNA",
                            id = "donor_id",
                            pred = c("ADNC", "sex", "PMI", "assay", "Ageatdeath", "APOE4_status","self_reported_ethnicity"),
                            offset = "nCount_RNA")

order_index <- order(neb_data$id)
# Reorder the count matrix by the id
neb_data$count <- neb_data$count[, order_index]
# Reorder the other components of neb_data
neb_data$id <- neb_data$id[order_index]
neb_data$pred <- neb_data$pred[order_index, ]
neb_data$offset <- neb_data$offset[order_index]

df <- model.matrix( ~ ADNC + sex + PMI + assay + Ageatdeath + APOE4_status + self_reported_ethnicity,
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
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData",
              "Applying NEBULA on SEA-AD.")

save(nebula_res, 
     date_of_run, session_info, note,
     start_time, end_time,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup6/Writeup6_SEA-AD_NEBULA.RData")

print("Done! :)")
