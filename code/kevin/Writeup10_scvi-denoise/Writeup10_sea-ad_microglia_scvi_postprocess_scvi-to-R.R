rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(arrow)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_cleaned.RData")
seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)

# Reading the Feather file
df <- arrow::read_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scvi.feather")
rowname_vec <- as.matrix(df[,ncol(df)])[,1]
df <- df[,-ncol(df)]
df <- as.matrix(df)
rownames(df) <- rowname_vec

# df2 <- scale(df); cor_mat <- crossprod(df2)/nrow(df2); diag(cor_mat) <- NA

Seurat::DefaultAssay(seurat_obj) <- "RNA"

Seurat::VariableFeatures(seurat_obj) <- colnames(df)
seurat_obj <- subset(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
SeuratObject::LayerData(object = seurat_obj, 
                        assay = "RNA", 
                        layer = "data") <- t(df)

# compute PCA and UMAP
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             features = Seurat::VariableFeatures(seurat_obj),
                             verbose = FALSE)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:30)

##############################

# fixing all the defficiencies
# update the RNA assay
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

# put factors
factor_vars_list <- list(Braakstage = TRUE,
                         CERADscore = c("Absent", "Sparse", "Frequent", "Moderate"),
                         LATENCstage = TRUE,
                         Microinfarctpathology = TRUE,
                         Thalphase = TRUE,
                         sex = TRUE,
                         assay = TRUE,
                         self_reported_ethnicity = TRUE,
                         development_stage = TRUE,
                         Cognitivestatus = c("Nodementia", "Dementia"),
                         ADNC = c("NotAD", "Low", "Intermediate", "High"),
                         APOE4status = c("N", "Y"),
                         donor_id = names(sort(table(seurat_obj$donor_id), decreasing = TRUE)),
                         seurat_clusters = TRUE,
                         integrated_snn_res.0.3 = TRUE,
                         genotype_APOE = TRUE,
                         orig.ident = TRUE,
                         Yearsofeducation = TRUE,
                         Supertype = TRUE,
                         Lewybodydiseasepathology = TRUE,
                         assay_ontology_term_id = TRUE,
                         suspension_type = TRUE,
                         cell_type_ontology_term_id = TRUE,
                         development_stage_ontology_term_id = TRUE,
                         disease_ontology_term_id = TRUE,
                         self_reported_ethnicity_ontology_term_id = TRUE,
                         organism_ontology_term_id = TRUE,
                         sex_ontology_term_id = TRUE,
                         tissue_ontology_term_id = TRUE,
                         is_primary_data = TRUE,
                         Neurotypicalreference = TRUE,
                         Class = TRUE,
                         Subclass = TRUE,
                         tissue_type = TRUE,
                         cell_type = TRUE,
                         disease = TRUE,
                         organism = TRUE,
                         tissue = TRUE)
for(kk in 1:length(factor_vars_list)){
  variable <- names(factor_vars_list)[kk]
  col_idx <- which(colnames(seurat_obj@meta.data) == variable)
  vec <- as.character(seurat_obj@meta.data[,col_idx])
  if(all(factor_vars_list[[kk]] == TRUE)){
    level_vec <- sort(unique(vec))
  } else {
    level_vec <- factor_vars_list[[kk]]
  }
  seurat_obj@meta.data[,col_idx] <- factor(vec, levels = level_vec)
}

#######################

var_vec <- c(donor_id = FALSE,
             Cognitivestatus = TRUE,
             ADNC = TRUE,
             Supertype = TRUE,
             sex = TRUE,
             self_reported_ethnicity = FALSE,
             assay = TRUE,
             APOE4status = TRUE,
             Thalphase = TRUE,
             CERADscore = TRUE,
             Braakstage = TRUE)

pdf("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup10/Writeup10_sea-ad_microglia_scvi_umap-covariates.pdf",
    onefile = T, width = 5, height = 5)

for(kk in 1:length(var_vec)){
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = "umap",
                           group.by = names(var_vec)[kk],
                           raster = TRUE)
  if(!var_vec[kk]) plot1 <- plot1 + Seurat::NoLegend()
  print(plot1)
}

dev.off()

#############

SeuratObject::LayerData(object = seurat_obj, 
                        assay = "RNA", 
                        layer = "scale.data") <- NULL

# switch gene names
RNA <- seurat_obj[["RNA"]] 
ensg_vec <- rownames(RNA@features@.Data)
mapping <- seurat_obj[["RNA"]]@meta.data
rownames(mapping) <- mapping$var.features
gene_vec <- mapping[ensg_vec, "feature_name"]
gene_vec <- as.character(gene_vec)
stopifnot(!any(duplicated(gene_vec)))
rownames(RNA@features@.Data) <- gene_vec
seurat_obj[["RNA"]] <- RNA
Seurat::VariableFeatures(seurat_obj) <- gene_vec

head(SeuratObject::Features(seurat_obj))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- "After running scVI."
save(seurat_obj,
     date_of_run, session_info, note,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData")


