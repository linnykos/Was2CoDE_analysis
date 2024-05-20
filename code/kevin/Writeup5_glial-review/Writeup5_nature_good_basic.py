# scVI environment

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import scanpy as sc
import scvi
import seaborn as sns
import torch
import pandas as pd
import gc
import os

# Set global settings
scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

print("Last run with scvi-tools version:", scvi.__version__)

file_path = "/home/users/kzlin/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.h5ad"
adata = ad.read_h5ad(file_path)

donor_ids_to_keep = ["D:10", "D:13", "D:16", "D:1", "D:20", # control (M, F, F, F, M). First 3 on Batch 2, last 2 on Batch 1
                     "D:11", "D:2", "D:3", "D:12", "D:4"] # cases (M, F, F, M, F). First 3 on Batch 1, last 2 on Batch 2

# Filter the `obs` DataFrame to include only these donors
filtered_obs = adata.obs[adata.obs["Pt_ID"].isin(donor_ids_to_keep)]

# Create a new AnnData object with the filtered cells
filtered_adata = adata[filtered_obs.index].copy()
filtered_adata.layers["counts"] = filtered_adata.X.copy()

pd.crosstab(filtered_adata.obs['SeqBatch'], filtered_adata.obs['CognitiveStatus'])

filtered_adata.obs["SeqBatch"] = filtered_adata.obs["SeqBatch"].astype('category')
filtered_adata.obs["Pt_ID"] = filtered_adata.obs["Pt_ID"].astype('category')
filtered_adata.obs["Sex"] = filtered_adata.obs["Sex"].astype('category')
filtered_adata.obs["Race"] = filtered_adata.obs["Race"].astype('category')

# Normalizing to median total counts
sc.pp.normalize_total(filtered_adata)
# Logarithmize the data
sc.pp.log1p(filtered_adata)

sc.pp.highly_variable_genes(
    filtered_adata,
    n_top_genes=5000,
    batch_key="Pt_ID",
    subset=True
)

filtered_adata

sc.tl.pca(filtered_adata)
sc.pp.neighbors(filtered_adata)
sc.tl.umap(filtered_adata, min_dist=0.3)

# Path to the directory where the plot will be saved
save_dir = os.path.expanduser("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5")
save_path = os.path.join(save_dir, "Writeup5_nature_good_basic_Pt_ID.png")

# Create the directory if it does not exist
os.makedirs(save_dir, exist_ok=True)

# Create the UMAP plot without displaying it
sc.pl.umap(
    filtered_adata,
    color=["Pt_ID"],
    frameon=False,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

save_path = os.path.join(save_dir, "Writeup5_nature_good_basic_SeqBatch.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    filtered_adata,
    color=["SeqBatch"],
    frameon=False,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

save_path = os.path.join(save_dir, "Writeup5_nature_good_basic_CognitiveStatus.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    filtered_adata,
    color=["CognitiveStatus"],
    frameon=False,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

print("Done! :)")