rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
set.seed(10)

print("Starting!")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup11/Writeup11_rosmap_scVI-postprocessed.RData") 
# head(seurat_obj@meta.data)

# construct the age vector
# age_vec <- as.numeric(seurat_obj$age_death)
# seurat_obj$age_death <- age_vec
# seurat_obj$APOEe4_status <- seurat_obj@meta.data[,"APOEe4_status"]
#colnames(seurat_obj@meta.data)

# gene_vec <- Seurat::VariableFeatures(seurat_obj@meta.data[["RNA"]])

# tmp <- paste0("ID_", as.character(seurat_obj@meta.data$Pt_ID))
# seurat_obj$Pt_ID <- factor(tmp)

# Reclassify ADpath levels
# seurat_obj$ADpath <- factor(ifelse(seurat_obj$ADpath == "yes", "Case", "Control"))
# summary(seurat_obj$ADpath)

# table(seurat_obj$Pt_ID, seurat_obj$ADpath)
# summary(seurat_obj$APOEe4_status)
# table(seurat_obj$Pt_ID, seurat_obj$APOEe4_status)

categorical_vars <- c("ADpath", "sex", "race", "batch", "APOEe4_status")
numerical_vars <- c("age_death", "pmi")
# 
# for (variable in categorical_vars) {
#   seurat_obj@meta.data[, variable] <- factor(seurat_obj@meta.data[, variable])
# }
# 
# for (variable in numerical_vars) {
#   seurat_obj@meta.data[, variable] <- as.numeric(as.character(seurat_obj@meta.data[, variable]))
# }

zz <- seurat_obj@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)
###########
# gene_vec <- Seurat::VariableFeatures(seurat_obj[["RNA"]])
# seurat_obj <- subset(seurat_obj, features = gene_vec)

##########################

neb_data <- nebula::scToNeb(obj = seurat_obj,
                            assay = "RNA",
                            id = "Pt_ID",
                            pred = c("ADpath", "sex", "pmi", "batch", "age_death", "APOEe4_status","race"),
                            offset = "nCount_RNA")

order_index <- order(neb_data$id)
# Reorder the count matrix by the id
neb_data$count <- neb_data$count[, order_index]
# Reorder the other components of neb_data
neb_data$id <- neb_data$id[order_index]
neb_data$pred <- neb_data$pred[order_index, ]
neb_data$offset <- neb_data$offset[order_index]

df <- model.matrix( ~ ADpath + sex + pmi + batch + age_death + APOEe4_status + race,
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
