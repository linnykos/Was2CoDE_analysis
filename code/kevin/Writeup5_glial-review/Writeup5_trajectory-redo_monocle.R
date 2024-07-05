# this version is to be run on Klone (HYAK)

rm(list=ls())
library(monocle3, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(SeuratWrappers)

# following https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

load("/gscratch/jayadevlab/projects/single_cell/brain/pu1_only/2021-03-31_22_samples/data/r_files/Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")

# head(ss_data_norm@meta.data)

set.seed(10)
Seurat::DefaultAssay(ss_data_norm) <- "integrated"
ss_data_norm <- Seurat::ScaleData(ss_data_norm, 
                                  features = Seurat::VariableFeatures(ss_data_norm))
ss_data_norm <- Seurat::RunPCA(ss_data_norm, 
                               features = Seurat::VariableFeatures(ss_data_norm),
                               verbose = FALSE)

set.seed(10)
ss_data_norm <- Seurat::RunUMAP(ss_data_norm, 
                                dims = 1:30)

ss_data_norm <- subset(ss_data_norm, 
                       integrated_snn_res.0.3 != "10")

Seurat::DimPlot(ss_data_norm,
                reduction = "umap",
                group.by = "integrated_snn_res.0.3")

tmp_seurat <- Seurat::CreateSeuratObject(counts = t(ss_data_norm[["pca"]]@cell.embeddings[,1:50]))

cds <- SeuratWrappers::as.cell_data_set(ss_data_norm)

set.seed(10)
cds <- monocle3::reduce_dimension(cds)

# hot-wire the UMAP in
SingleCellExperiment::reducedDim(cds, type = "UMAP") <- ss_data_norm[["umap"]]@cell.embeddings

set.seed(10)
k_value <- 58
cds <- monocle3::cluster_cells(cds, k = k_value)
length(unique(cds@clusters[["UMAP"]]$clusters))

set.seed(10)
cds <- monocle3::learn_graph(cds,
                             use_partition = FALSE)

# hot-wire the clustering in
tmp <- ss_data_norm@meta.data[,"integrated_snn_res.0.3"]
tmp <- droplevels(tmp)
tmp <- as.numeric(as.character(tmp))
tmp2 <- tmp
# swap: current 4->new 6, current 6->new 7, current 7->4
tmp2[tmp == 4] <- 6
tmp2[tmp == 6] <- 7
tmp2[tmp == 7] <- 4
tmp2 <- factor(tmp2)
names(tmp2) <- Seurat::Cells(ss_data_norm)
cds@clusters[["UMAP"]]$clusters <- tmp2

monocle3::plot_cells(cds, 
                     label_groups_by_cluster = FALSE, 
                     label_leaves = FALSE, 
                     label_branch_points = FALSE)

#### hot-wire the thing

set.seed(10)
cds <- monocle3::learn_graph(cds,
                             use_partition = FALSE)

monocle3::plot_cells(cds, 
                     label_groups_by_cluster = FALSE, 
                     label_leaves = FALSE, 
                     label_branch_points = FALSE)

