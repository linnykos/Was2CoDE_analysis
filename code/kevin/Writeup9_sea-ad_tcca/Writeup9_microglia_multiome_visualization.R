rm(list=ls())

library(Seurat)
library(Signac)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup9/Writeup9_SEAAD_MTG_Microglia-PVM_multiome.RData")

categorical_vec <- c("Cognitive status",
                     "APOE4 status",
                     "LATE-NC stage", "CERAD score")
numerical_vec <- c("Continuous Pseudo-progression Score",
                   "Age at death", "Braak stage",
                   "Thal phase")

# fix covariates
variable <- "Braak stage"
tmp <- seurat_obj@meta.data[,variable]
tmp2 <- sapply(as.character(tmp), function(x){
  strsplit(x, split = " ")[[1]][2]
})
tmp2[tmp2 == "I"] <- "1"; tmp2[tmp2 == "II"] <- "2"
tmp2[tmp2 == "III"] <- "3"; tmp2[tmp2 == "IV"] <- "4"
tmp2[tmp2 == "V"] <- "5"; tmp2[tmp2 == "VI"] <- "6"
seurat_obj@meta.data[,variable] <- as.numeric(tmp2)

variable <- "Thal phase"
tmp <- seurat_obj@meta.data[,variable]
tmp2 <- sapply(as.character(tmp), function(x){
  strsplit(x, split = " ")[[1]][2]
})
seurat_obj@meta.data[,variable] <- as.numeric(tmp2)

variable <- "LATE-NC stage"
tmp <- seurat_obj@meta.data[,variable]
levels(tmp)[which(levels(tmp) == "Staging Precluded by FTLD with TDP43 or ALS/MND or TDP-43 pathology is unclassifiable")] <- "Unclassifiable"
seurat_obj@meta.data[,variable] <- tmp

variable <- "Age at death"
tmp <- seurat_obj@meta.data[,"development_stage"]
tmp <- as.character(tmp)
tmp[which(tmp == "80 year-old and over human stage")] <- "90-year-old human stage"
tmp2 <- sapply(tmp, function(x){
  strsplit(x, split = "-")[[1]][1]
})
seurat_obj@meta.data[,variable] <- as.numeric(tmp2)

# check covariates
for(variable in categorical_vec){
  print(variable)
  print(length(which(is.na(seurat_obj@meta.data[,variable]))))
  print(table(seurat_obj@meta.data[,variable]))
  print("====")
}

for(variable in numerical_vec){
  print(variable)
  print(length(which(is.na(seurat_obj@meta.data[,variable]))))
  print(quantile(seurat_obj@meta.data[,variable]))
  print("====")
}

seurat_obj[["umap_scVI"]] <- seurat_obj[["umap"]]
seurat_obj[["umap_MultiVI"]] <- seurat_obj[["umap_atac"]]

#################

# create color gradient for donors
tab_mat <- table(seurat_obj$donor_id, seurat_obj@meta.data[,"Cognitive status"])
tab_mat <- tab_mat[,c("No dementia", "Dementia")]
tab_mat <- tab_mat[rowSums(tab_mat)!=0,]

ad_donor <- rownames(tab_mat)[which(tab_mat[,"Dementia"] != 0)]
ctrl_donor <- rownames(tab_mat)[which(tab_mat[,"Dementia"] == 0)]

case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(length(ad_donor))
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(length(ctrl_donor))
col_vec <- c(case_color_palette, control_color_palette)
names(col_vec) <- c(ad_donor, ctrl_donor)

#################

# make some clusterings

set.seed(10)
Seurat::DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, 
                                           selection.method = "vst", nfeatures = 2000)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, 
                             features = Seurat::VariableFeatures(seurat_obj),
                             verbose = FALSE)
set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)

set.seed(10)
Seurat::DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- Signac::RunTFIDF(seurat_obj)
seurat_obj <- Signac::FindTopFeatures(seurat_obj, min.cutoff = "q0")
seurat_obj <- Signac::RunSVD(seurat_obj)
set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, 
                                    reduction = 'lsi',
                                    dims = 2:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                              reduction = "lsi", 
                              dims = 2:30, 
                              reduction.name = "umap_atac", 
                              reduction.key = "atacUMAP_")

seurat_obj <- Seurat::FindMultiModalNeighbors(
  seurat_obj, 
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:30), 
  modality.weight.name = "RNA.weight"
)
seurat_obj <- Seurat::RunUMAP(seurat_obj, nn.name = "weighted.nn", 
                              reduction.name = "wnn.umap", 
                              reduction.key = "wnnUMAP_")

# create plots
plot_filename <- c("RNA", "ATAC", "WNN", "scVI", "MultiVI")
umap_name <- c("umap", "umap_atac", "wnn.umap" ,"umap_scVI", "umap_MultiVI")

for(kk in 1:length(plot_filename)){
  print(plot_filename[kk])
  
  file <- plot_filename[kk]
  umap <- umap_name[kk]
  
  pdf(paste0("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup9/Writeup9_microglia_", file,"_plots.pdf"),
      onefile = T, width = 7, height = 5)
  
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = umap,
                           group.by = "donor_id", 
                           cols = col_vec)
  print(plot1)
  
  for(variable in c(categorical_vec, "RNA_snn_res.0.5", "ATAC_snn_res.0.5")){
    plot1 <- Seurat::DimPlot(seurat_obj, 
                             reduction = umap,
                             group.by = variable)
    print(plot1)
  }
  
  for(variable in numerical_vec){
    plot1 <- Seurat::FeaturePlot(seurat_obj, 
                                 reduction = umap,
                                 features = variable)
    print(plot1)
  }
  
  dev.off()
}
