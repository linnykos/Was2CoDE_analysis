rm(list=ls())
library("Seurat")

seurat_obj <- readRDS("~/kzlinlab/data/sea-ad/microglia-pvm_dpc.rds")

options(Seurat.object.assay.version = "v5")
seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

#########

seurat_obj$CognitiveStatus <- seurat_obj$`Cognitive status`  
seurat_obj$`Cognitive status` <- NULL
seurat_obj$AgeAtDeath <- seurat_obj$`Age at death`  
seurat_obj$`Age at death` <- NULL
seurat_obj$YearsOfEducation <- seurat_obj$`Years of education`  
seurat_obj$`Years of education` <- NULL
seurat_obj$FractionMitochrondrialUMIs <- seurat_obj$`Fraction mitochrondrial UMIs`  
seurat_obj$`Fraction mitochrondrial UMIs` <- NULL

# remove all spaces in all covariates
meta_vars <- colnames(seurat_obj@meta.data)
for(variable in meta_vars){
  vec <- seurat_obj@meta.data[,variable]
  if(is.character(vec) | is.factor(vec)){
    vec <- as.character(vec)
    vec <- gsub(pattern = "[^[:alnum:] ]", replacement = "", x = vec)
    vec <- gsub(pattern = " ", replacement = "", x = vec)
    seurat_obj@meta.data[,variable] <- factor(vec)
  }
}

############

rownames(seurat_obj[["RNA"]]) <- seurat_obj[["RNA"]]@meta.data[,"feature_name"]

Seurat::DefaultAssay(seurat_obj) <- "RNA"
set.seed(10)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                           selection.method = "vst",
                                           nfeatures = 2000)

sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
) 
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
sun_genes <- intersect(sun_genes, SeuratObject::Features(seurat_obj))

prater_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/43587_2023_424_MOESM3_ESM.xlsx",
  sheet = "DEGs_cluster_to_all",
  startRow = 2
) 
prater_list <- list(
  vs_all = intersect(
    unique(prater_sheet$gene[which(prater_sheet$p_val_adj <= 0.05)]), 
    SeuratObject::Features(seurat_obj))
)
prater_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/43587_2023_424_MOESM3_ESM.xlsx",
  sheet = "DEGs_cluster_to_clust1",
  startRow = 2
) 
prater_list[["vs_1"]] <- intersect(
  unique(prater_sheet$Gene[which(prater_sheet$p_val_adj <= 0.05)]),
  SeuratObject::Features(seurat_obj)
)
prater_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/43587_2023_424_MOESM3_ESM.xlsx",
  sheet = "AD_vs_Ctrl_DEGs",
  startRow = 3
) 
prater_list[["AD_vs_Ctrl"]] <- intersect(
  unique(prater_sheet$Gene[which(prater_sheet$padj <= 0.05)]),
  SeuratObject::Features(seurat_obj)
)

#################

hk_df <- read.csv("~/kzlinlab/projects/subject-de/data/HSIAO_HOUSEKEEPING_GENES.v2023.2.Hs.tsv", 
                   sep = "\t")
msigdb_hk <- hk_df[which(hk_df$STANDARD_NAME == "GENE_SYMBOLS"),"HSIAO_HOUSEKEEPING_GENES"]
msigdb_hk <- strsplit(msigdb_hk, split = ",")[[1]]
msigdb_hk <- msigdb_hk[which(sapply(msigdb_hk, nchar) > 0)]
msigdb_hk <- intersect(
  msigdb_hk,
  SeuratObject::Features(seurat_obj)
)

hk_df <- read.csv("~/kzlinlab/projects/subject-de/data/Housekeeping_GenesHuman.csv",
                  sep = ";")
hounkpe_hk <- hk_df[,"Gene.name"]
hounkpe_hk <- intersect(
  hounkpe_hk,
  SeuratObject::Features(seurat_obj)
)

#################

seurat_obj@misc <- list(
  sun_genes = sun_genes,
  prater_list = prater_list,
  msigdb_hk = msigdb_hk,
  hounkpe_hk = hounkpe_hk
)

Seurat::VariableFeatures(seurat_obj) <- sort(unique(
  c(Seurat::VariableFeatures(seurat_obj), 
    sun_genes,
    unlist(prater_list),
    msigdb_hk, 
    hounkpe_hk)
))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Working from ~/kzlinlab/data/sea-ad/microglia-pvm_dpc.rds.",
              "Modified the metadata names to not be weird, and added Variable Features,",
              "two sets from microglia papers, two sets from housekeeping papers")

save(date_of_run, session_info, note,
     seurat_obj,
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup2/Writeup2_sea-ad_microglia_preprocess.RData")

