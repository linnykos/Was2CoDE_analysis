rm(list=ls())
library(Seurat)

data_folder <- "~/kzlinlab/data/synapse_ad_rosmap-microglia/"
counts <- readRDS(paste0(data_folder, "ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds"))
metadata_cells <- readRDS(paste0(data_folder, "ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds"))
donor_codebook <- read.csv(paste0(data_folder, "MIT_ROSMAP_Multiomics_individual_metadata.csv"))
metadata_donor <- read.csv(paste0(data_folder, "ROSMAP_clinical.csv"))
celltype_seurat_obj <- readRDS(paste0(data_folder, "Immune_cells.rds"))
metadata_donor2 <- read.csv(paste0(data_folder, "individual_metadata_deidentified.tsv"), 
                            sep = "\t")

#######

# build off of metadata_cells
# First: transfer over celltype labels
barcode_vec <- sapply(rownames(metadata_cells), function(x){
  tmp <- strsplit(x, split = "\\.")[[1]][2]
  tmp <- strsplit(tmp, split = "-")[[1]][1]
  tmp
})
celltype_seurat_barcodes <- Seurat::Cells(celltype_seurat_obj)
celltype_seurat_barcodes <- sapply(celltype_seurat_barcodes, function(x){
  strsplit(x, split = "-")[[1]][1]
})
celltype_vec <- celltype_seurat_obj$cell_type_high_resolution
names(celltype_vec) <- celltype_seurat_barcodes
celltype_vec_transferred <- celltype_vec[barcode_vec]
celltype_vec_transferred <- as.character(celltype_vec_transferred)
celltype_vec_transferred <- gsub(pattern = " ", 
                                 replace = "_",
                                 celltype_vec_transferred)
celltype_vec_transferred <- factor(celltype_vec_transferred)
metadata_cells$celltype <- celltype_vec_transferred

# Enumerate all the donor's ROSMAP ID and individualID
donor_codebook <- donor_codebook[which(!duplicated(donor_codebook$subject)),]
rownames(donor_codebook) <- donor_codebook$subject
donor_names <- sort(unique(metadata_cells$subject))
individual_names <- donor_codebook[donor_names, "individualID"]
names(individual_names) <- donor_names
metadata_cells$Pt_ID <- individual_names[metadata_cells$subject]

# do some simple fixing
metadata_donor$race[which(metadata_donor$race == "(Other)")] <- 6
metadata_donor$race[is.na(metadata_donor$race)] <- 7
vec <- which(!is.na(metadata_donor$apoe_genotype))
metadata_donor$apoe_genotype[vec] <- paste0("apoe", metadata_donor$apoe_genotype[vec])

# Convert relevant covariates in metadata_donor into factors
categorical_vars <- c("Study", "msex", "race", "spanish", "apoe_genotype", "projid", "individualID")
for(variable in categorical_vars){
  metadata_donor[,variable] <- factor(metadata_donor[,variable])
}
# summary(metadata_donor)

age_variables <- c("age_at_visit_max", "age_first_ad_dx", "age_death")
for(age_variable in age_variables){
  vec <- metadata_donor[,age_variable] 
  vec[which(vec == "90+")] <- "90"
  metadata_donor[,age_variable] <- as.numeric(vec)
}

msex_levels <- c(Female = 0, Male = 1)
race_levels <- c(White = 1, 
                 BlackAfrican = 2, 
                 IndianNative= 3, 
                 HawaiianPacific = 4, 
                 Asian = 5,
                 Other = 6,
                 Unknown = 7)
spanish_levels <- c(YesOrigins = 1, NoOrigins = 2)

metadata_donor$msex <- factor(names(msex_levels)[metadata_donor$msex])
metadata_donor$race <- factor(names(race_levels)[metadata_donor$race])
metadata_donor$spanish <- factor(names(spanish_levels)[metadata_donor$spanish])
metadata_donor$sex <- metadata_donor$msex
metadata_donor$APOEe4_status <- factor(as.character(metadata_donor$apoe_genotype %in% c("apoe34", "apoe44")))
vec <- as.character(metadata_donor$Study)
vec[vec == ""] <- NA
metadata_donor$Study <- factor(vec)
rownames(metadata_donor) <- metadata_donor$individualID

# start transferring variables over
categorical_variables <- c("Study", "race", "spanish", "apoe_genotype", 
                           "sex", "APOEe4_status")
numerical_variables <- c("educ", "age_at_visit_max", "age_first_ad_dx",
                         "age_death", "cts_mmse30_first_ad_dx", 
                         "cts_mmse30_lv", "pmi", 
                         "braaksc", "ceradsc", "cogdx", "dcfdx_lv")
for(variable in categorical_variables){
  metadata_cells[,variable] <- factor(metadata_donor[metadata_cells$Pt_ID, variable])
}
for(variable in numerical_variables){
  metadata_cells[,variable] <- metadata_donor[metadata_cells$Pt_ID, variable]
}

# do some extra formatting
metadata_cells$msex <- NULL
categorical_variables <- c("subject", "brainRegion", "batch", 
                           "ADdiag3types", "Pt_ID")
for(variable in categorical_variables){
  metadata_cells[,variable] <- factor(metadata_cells[,variable])
}
metadata_cells <- metadata_cells[colnames(counts),]

# lastly, add the AD pathology
rownames(metadata_donor2) <- metadata_donor2$subject
metadata_cells$ADpath <- factor(metadata_donor2[metadata_cells$subject, "Pathologic_diagnosis_of_AD"],
                                levels = c("no", "yes"))
table(metadata_cells$ADdiag3types, metadata_cells$ADpath)

##############################

seurat_obj <- Seurat::CreateSeuratObject(counts = counts,
                                         meta.data = metadata_cells)
# remove the subject with no pathology
keep_vec <- rep(TRUE, length(Seurat::Cells(seurat_obj)))
keep_vec[is.na(seurat_obj$ADpath)] <- FALSE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# remove all non-microglia
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[grep("Mic", seurat_obj$celltype)] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# keep only PFC
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[which(seurat_obj$brainRegion == "PFC")] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# rename any categoricals to be not numbers
seurat_obj$seurat_clusters <- factor(paste0("c", seurat_obj$seurat_clusters))

# remove any person with less than 50 cells
donor_cellnumbers <- table(seurat_obj$subject)
donor_keep <- names(donor_cellnumbers)[which(donor_cellnumbers >= 50)]
keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_obj)))
keep_vec[seurat_obj$subject %in% donor_keep] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

categorical_variables <- c("subject", "brainRegion", "batch", 
                           "ADdiag3types", "Pt_ID", "Study", 
                            "race", "spanish", "apoe_genotype", 
                            "sex", "APOEe4_status", "celltype",
                           "seurat_clusters")
for(variable in categorical_variables){
  seurat_obj@meta.data[,variable] <- droplevels(seurat_obj@meta.data[,variable])
}

seurat_obj
summary(seurat_obj@meta.data)
SeuratObject::LayerData(seurat_obj,
                        assay = "RNA", 
                        layer = "counts")[1:10,1:10]

####

date_of_run <- Sys.time()
session_info <- devtools::session_info()

print("Saving")
out_folder <- "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup11/"
save(seurat_obj, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11_rosmap-initial.RData"))






