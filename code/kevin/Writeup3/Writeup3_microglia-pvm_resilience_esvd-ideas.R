rm(list=ls())
library(Seurat)
library(ideas)
library(foreach)
library(doRNG)
RNGkind("L'Ecuyer-CMRG")
library(doParallel)
doParallel::registerDoParallel(cores=6)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_esvd.RData")
set.seed(10)

# process the seurat_obj$donor_id
donor_vec <- seurat_obj$donor_id
donor_vec <- sapply(donor_vec, function(x){
  paste0(substr(x, start = 0, stop = 3), ".", 
         substr(x, start = 4, stop = 5), ".",
         substr(x, start = 6, stop = 8))
})
seurat_obj$donor_id <- donor_vec

# remove the reference cells
keep_vec <- rep(FALSE, length(SeuratObject::Cells(seurat_obj)))
names(keep_vec) <- SeuratObject::Cells(seurat_obj)
keep_vec[rownames(eSVD_obj$dat)] <- TRUE
seurat_obj$keep <- keep_vec
seurat_obj <- subset(seurat_obj, keep == TRUE)

# add the resiliency phenotype back in
stopifnot(all(SeuratObject::Cells(seurat_obj) == rownames(eSVD_obj$dat)))
resiliency_vec <- eSVD_obj$covariates[,"resilient_Resilient"]
seurat_obj$resiliency <- resiliency_vec

#########################

# see example code:
# https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1d_saver.R#L177 on how rd was defined
# https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1c_ideas.R for how the method was used

# notes:
## the only things in meta_cell that are needed is: "cell_id", "individual", and var_per_cell (which is the read depth)
## meta_ind must have a column called "individual" and var2test
## colnames in mat must be equal to meta_cell$cell_id 
## see: https://github.com/Sun-lab/ideas/blob/main/R/ideas_dist.R

cell_covariates <- data.frame(cell_id = rownames(seurat_obj@meta.data),
                              individual = seurat_obj$donor_id,
                              rd = seurat_obj$nCount_RNA)
tab_mat <- table(seurat_obj$donor_id, seurat_obj$resiliency)
bool_vec <- rep(FALSE, nrow(tab_mat))
bool_vec[which(tab_mat[,"1"] > 0)] <- TRUE
indiv_covariates <- data.frame(individual = rownames(tab_mat),
                               resiliency = bool_vec)

#######################

# ideas::ideas_dist can only take in dense matrices, so we're going to limit
#   which genes we analyze for now

prater_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/43587_2023_424_MOESM3_ESM.xlsx",
  sheet = "AD_vs_Ctrl_DEGs",
  startRow = 3
) 
prater_genes <- sort(unique(prater_sheet[which(prater_sheet[,"padj"] <= 0.05),"Gene"]))
prater_genes <- intersect(prater_genes, Seurat::VariableFeatures(seurat_obj))

mat <- SeuratObject::LayerData(seurat_obj, 
                               assay = "RNA", 
                               features = prater_genes,
                               layer = "counts")
mat <- as.matrix(mat)
stopifnot(all(colnames(mat) == cell_covariates$cell_id))

print("Starting ideas::ideas_dist")
set.seed(10)
time_start <- Sys.time()
dist_tensor <- ideas::ideas_dist(count_input = mat,
                                 meta_cell = cell_covariates,
                                 meta_ind = indiv_covariates,
                                 var_per_cell = "rd",
                                 var2test = "resiliency",
                                 var2test_type = "binary",
                                 d_metric = "Was",
                                 fit_method = "nb")
time_end <- Sys.time()

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(date_of_run, session_info, 
     time_start, time_end,
     cell_covariates, dist_tensor, indiv_covariates,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_ideas.RData")

print("Done! :)")