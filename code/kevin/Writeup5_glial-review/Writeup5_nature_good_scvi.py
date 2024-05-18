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

scvi.model.SCVI.setup_anndata(
    filtered_adata,
    layer="counts",
    categorical_covariate_keys=["Sex", "Pt_ID", "Race"],
    continuous_covariate_keys=["percent.mito"],
    batch_key="SeqBatch"
)

model = scvi.model.SCVI(filtered_adata,
                        n_layers=2,
                        n_latent=30,
                        gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
filtered_adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

model.save(dir_path="/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_scvi-model", 
           overwrite=True)

reserved_names = {'_index'}
if filtered_adata.raw is not None:
    raw_var_reserved = reserved_names.intersection(filtered_adata.raw.var.columns)
    print("Reserved names in raw.var:", raw_var_reserved)
    if '_index' in filtered_adata.raw.var.columns:
        filtered_adata.raw.var.rename(columns={'_index': 'index_raw_var'}, inplace=True)

filtered_adata.write("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup5/Writeup5_nature_good_anndata.h5ad")
