# scVI environment

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import scanpy as sc
import scvi
import seaborn as sns
import torch
import gc

# Set global settings
scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

print("Last run with scvi-tools version:", scvi.__version__)

file_path = "/home/users/kzlin/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.h5ad"
adata = ad.read_h5ad(file_path)

donor_ids_to_keep = ["D:14", "D:20", "D:21", "D:18", "D:22", "D:13", "D:4", "D:1", "D:12", "D:11"] 

# Filter the `obs` DataFrame to include only these donors
filtered_obs = adata.obs[adata.obs["Pt_ID"].isin(donor_ids_to_keep)]

# Create a new AnnData object with the filtered cells
filtered_adata = adata[filtered_obs.index].copy()

filtered_adata.obs["SeqBatch"] = filtered_adata.obs["SeqBatch"].astype('category')
filtered_adata.obs["Pt_ID"] = filtered_adata.obs["Pt_ID"].astype('category')
filtered_adata.obs["Sex"] = filtered_adata.obs["Sex"].astype('category')

sc.pp.highly_variable_genes(
    filtered_adata,
    n_top_genes=2000,
    batch_key="Pt_ID",
    subset=True
)

filtered_adata

scvi.model.SCVI.setup_anndata(
    filtered_adata,
    categorical_covariate_keys=["Sex", "Pt_ID"],
    continuous_covariate_keys=["percent.mito"],
    batch_key="SeqBatch"
)

model = scvi.model.SCVI(filtered_adata)
model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

model.save(dir_path="/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_scvi-model", 
           overwrite=True)
