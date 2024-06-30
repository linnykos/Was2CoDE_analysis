# scVI environment

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import scanpy as sc
import numpy as np
import scvi
import seaborn as sns
import torch
import gc
import os
import sys
import pandas as pd

# Set global settings
scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

print("Last run with scvi-tools version:", scvi.__version__)

adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-noCovariates_anndata.h5ad")

# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
# use scVI latent space for UMAP generation
SCVI_LATENT_KEY = "X_scVI"
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)

# Path to the directory where the plot will be saved
save_dir = os.path.expanduser("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup6")
# Create the directory if it does not exist
os.makedirs(save_dir, exist_ok=True)

# Generate a shuffled index array
shuffled_indices = np.random.permutation(adata.n_obs)
# Reindex the AnnData object
adata2 = adata[shuffled_indices, :]

save_path = os.path.join(save_dir, "Writeup6_prater_scvi-noCovariates_Pt_ID.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    adata2,
    color=["Pt_ID"],
    frameon=False,
    title="By Donor",
    size=5,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

save_path = os.path.join(save_dir, "Writeup6_prater_scvi-noCovariates_CognitiveStatus.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    adata2,
    color=["CognitiveStatus"],
    frameon=False,
    title="By Cognitive Status",
    size=5,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

save_path = os.path.join(save_dir, "Writeup6_prater_scvi-noCovariates_SeqBatch.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    adata2,
    color=["SeqBatch"],
    frameon=False,
    title="By Batch",
    size=5,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

save_path = os.path.join(save_dir, "Writeup6_prater_scvi-noCovariates_SeuratClusters.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    adata2,
    color=["seurat_clusters"],
    frameon=False,
    title="By Seurat Clusters",
    size=5,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

##############

adata.obs["Pt_ID"] = adata.obs["Pt_ID"].astype('category')

# denoise the gene expression
# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
model = scvi.model.SCVI.load("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-noCovariates-model", 
                             adata)

unique_batches = adata.obs["Pt_ID"].unique().tolist()
denoised = model.get_normalized_expression(adata, 
                                           return_mean=True, 
                                           transform_batch=unique_batches,
                                           library_size=1e4)

# Assuming df is your DataFrame
denoised.to_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-noCovariates-denoised.feather")