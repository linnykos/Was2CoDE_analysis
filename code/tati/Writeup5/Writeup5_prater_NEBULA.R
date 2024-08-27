rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")

gene_vec <- Seurat::VariableFeatures(ss_data_norm[["RNA"]])
ss_data_norm <- subset(ss_data_norm, features = gene_vec)

# adjust the covariates
# need to convert to numeric
tmp <- paste0("ID_", as.character(ss_data_norm$Pt_ID))
ss_data_norm$Pt_ID <- factor(tmp)

categorical_vars <- c("Sex", "SeqBatch", "Race", "Study_Designation","APOEe4_status")
numerical_vars <- c("PMI","coded_Age")

tmp <- as.character(ss_data_norm$coded_Age)
ss_data_norm$coded_Age <- tmp

for(variable in categorical_vars){
  ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
}
for(variable in numerical_vars){
  ss_data_norm@meta.data[,variable] <- as.numeric(as.character(ss_data_norm@meta.data[,variable]))
}

zz <- ss_data_norm@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)

##########################

neb_data <- nebula::scToNeb(obj = ss_data_norm,
                            assay = "RNA",
                            id = "Pt_ID",
                            pred = c("Study_Designation", "Sex", "PMI", "SeqBatch", "coded_Age", "APOEe4_status","Race"),
                            offset = "nCount_RNA")
df <- model.matrix( ~ Study_Designation + Sex + PMI + SeqBatch + coded_Age + APOEe4_status + Race,
                    data = neb_data$pred)
head(df)
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
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData.",
              "Applying NEBULA.")

save(nebula_res, 
     date_of_run, session_info, note,
     start_time, end_time,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_NEBULA.RData")

print("Done! :)")
