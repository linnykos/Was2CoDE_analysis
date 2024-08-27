rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
library(nebula)
set.seed(10)

print("Starting")

load("~/kzlinlab/projects/subject-de/out/kevin/preprocess/naive-preprocess.RData")

gene_vec <- Seurat::VariableFeatures(ss_data_norm[["RNA"]])
ss_data_norm <- subset(ss_data_norm, features = gene_vec)

# adjust the covariates
# need to convert to numeric
tmp <- paste0("ID_", as.character(ss_data_norm$Pt_ID))
ss_data_norm$Pt_ID <- factor(tmp)

categorical_vars <- c("Sex", "SeqBatch", "Race", "Study_Designation","genotype_APOE")
numerical_vars <- c("PMI","coded_Age")

tmp <- as.character(ss_data_norm$coded_Age)
tmp[which(tmp == "90+")] <- "90"
ss_data_norm$coded_Age <- tmp

for(variable in categorical_vars){
  ss_data_norm@meta.data[,variable] <- factor(ss_data_norm@meta.data[,variable])
}
for(variable in numerical_vars){
  ss_data_norm@meta.data[,variable] <- as.numeric(as.character(ss_data_norm@meta.data[,variable]))
}

# need to fill in NAs in PMI,Race
tab_list <- lapply(categorical_vars, function(variable){
  table(ss_data_norm$Pt_ID, ss_data_norm@meta.data[,variable])
})
names(tab_list) <- categorical_vars

pmi_missing_values <- c(ID_1 = 5.08,
                       ID_10 = 3.33,
                       ID_3 = 5.08,
                       ID_15 = 6.97)
for(i in 1:length(pmi_missing_values)){
  subj_id <- names(pmi_missing_values)[i]
  idx <- which(ss_data_norm$Pt_ID == subj_id)
  ss_data_norm$PMI[idx] <- pmi_missing_values[i]
}

race_missing_values <- c(ID_5 = "White")
for(i in 1:length(race_missing_values)){
  subj_id <- names(race_missing_values)[i]
  idx <- which(ss_data_norm$Pt_ID == subj_id)
  ss_data_norm$Race[idx] <- race_missing_values[i]
}

zz <- ss_data_norm@meta.data[,c(categorical_vars, numerical_vars)]
stopifnot(!any(is.na(zz)))
summary(zz)

##########################

neb_data <- nebula::scToNeb(obj = ss_data_norm,
                            assay = "RNA",
                            id = "Pt_ID",
                            pred = c("Study_Designation", "Sex", "PMI", "SeqBatch", "coded_Age", "genotype_APOE","Race"),
                            offset = "nCount_RNA")
df <- model.matrix( ~ Study_Designation + Sex + PMI + SeqBatch + coded_Age + genotype_APOE + Race,
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
note <- paste("Working from ~/kzlinlab/projects/subject-de/out/kevin/preprocess/naive-preprocess.RData.",
              "Applying NEBULA.")

save(nebula_res, 
     date_of_run, session_info, note,
     start_time, end_time,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup3/Writeup3_NEBULA.RData")

print("Done! :)")
