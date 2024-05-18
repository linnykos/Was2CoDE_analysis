# scVI environment

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import scanpy as sc
import scvi
import seaborn as sns
import torch
import gc
import os

# Set global settings
scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

print("Last run with scvi-tools version:", scvi.__version__)

filtered_adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_bad_anndata.h5ad")

# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
# use scVI latent space for UMAP generation
SCVI_LATENT_KEY = "X_scVI"
sc.pp.neighbors(filtered_adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(filtered_adata, min_dist=0.3)

sc.pl.umap(
    filtered_adata,
    color=["Pt_ID"],
    frameon=False,
    show=False
)

# Path to the directory where the plot will be saved
save_dir = os.path.expanduser("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup5")
save_path = os.path.join(save_dir, "Writeup5_nature_bad_scvi_Pt_ID.png")

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

save_path = os.path.join(save_dir, "Writeup5_nature_bad_scvi_SeqBatch.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    filtered_adata,
    color=["SeqBatch"],
    frameon=False,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')

save_path = os.path.join(save_dir, "Writeup5_nature_bad_scvi_CognitiveStatus.png")
# Create the UMAP plot without displaying it
sc.pl.umap(
    filtered_adata,
    color=["CognitiveStatus"],
    frameon=False,
    show=False
)
# Save the figure manually
plt.savefig(save_path, bbox_inches='tight')
