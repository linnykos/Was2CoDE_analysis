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

adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi2_anndata.h5ad")

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

save_path = os.path.join(save_dir, "Writeup6_prater_scvi2_Pt_ID.png")
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

save_path = os.path.join(save_dir, "Writeup6_prater_scvi2_CognitiveStatus.png")
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

save_path = os.path.join(save_dir, "Writeup6_prater_scvi2_SeqBatch.png")
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

save_path = os.path.join(save_dir, "Writeup6_prater_scvi2_SeuratClusters.png")
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

adata.obs["SeqBatch"] = adata.obs["SeqBatch"].astype('category')
adata.obs["Pt_ID"] = adata.obs["Pt_ID"].astype('category')
adata.obs["Sex"] = adata.obs["Sex"].astype('category')
adata.obs["Race"] = adata.obs["Race"].astype('category')

# denoise the gene expression
# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
model = scvi.model.SCVI.load("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi2-model", 
                             adata)

# following https://discourse.scverse.org/t/using-categorical-covariate-keys-when-sampling-or-generating-normalised-expression/1493
# manually set the mode and mean
categorical_covariates = ["Sex", "Race", "Pt_ID"]
# Store the original categories
original_categories = {covariate: adata.obs[covariate].cat.categories for covariate in categorical_covariates}

# Calculate mode for categorical covariates
for covariate in categorical_covariates:
    mode_value = adata.obs[covariate].mode()[0]
    adata.obs[covariate] = mode_value
    # Ensure the categorical variable type and categories are retained
    adata.obs[covariate] = adata.obs[covariate].astype('category')
    adata.obs[covariate] = adata.obs[covariate].cat.set_categories(original_categories[covariate])

# Calculate mean for continuous covariates
continuous_covariates = ["percent.mito", "coded_Age"]
for covariate in continuous_covariates:
    mean_value = adata.obs[covariate].mean()
    adata.obs[covariate] = mean_value

unique_batches = adata.obs["SeqBatch"].unique().tolist()
denoised = model.get_normalized_expression(adata, 
                                           return_mean=True, 
                                           transform_batch=unique_batches)

# Assuming df is your DataFrame
denoised.to_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-denoised.feather")