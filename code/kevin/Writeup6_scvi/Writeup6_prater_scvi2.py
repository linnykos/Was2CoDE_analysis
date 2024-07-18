# scVI environment

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import scanpy as sc
import scvi
import seaborn as sns
import torch
import pandas as pd
import numpy as np
import gc
from collections import Counter

# Set global settings
scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

print("Last run with scvi-tools version:", scvi.__version__)

file_path = "/home/users/kzlin/kzlinlab/data/microglia-prater-2023/Prater_Green_PU1_MGsubset_10clusters_DeID.h5ad"
adata = ad.read_h5ad(file_path)

# Counter(adata.obs["Sex"])
# Counter(adata.obs["coded_Age"])
# Counter(adata.obs["CognitiveStatus"])
# Counter(adata.obs["genotype_APOE"])
# Counter(adata.obs["Study_Designation"])
# adata.obs["PMI"].quantile([0, 0.25, 0.5, 0.75, 1])
# adata.obs["coded_Age"].quantile([0, 0.25, 0.5, 0.75, 1])
# adata.obs["BrainPh"].quantile([0, 0.25, 0.5, 0.75, 1])

adata.obs["SeqBatch"] = adata.obs["SeqBatch"].astype('category')
adata.obs["Sex"] = adata.obs["Sex"].astype('category')
adata.obs["Race"] = adata.obs["Race"].astype('category')
adata.obs["genotype_APOE"] = adata.obs["genotype_APOE"].astype('category')
adata.obs["CERAD"] = adata.obs["CERAD"].astype('category')

# adding raw counts for referring to it in the future
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    batch_key="Pt_ID",
    subset=True
)

adata

# # Check for NAs in the observation metadata (adata.obs)
# na_in_obs = adata.obs.isnull().values.any()
# print(f"Are there NAs in the observation metadata (adata.obs)? {na_in_obs}")

# # Identify columns in adata.obs with NA values
# na_columns_obs = adata.obs.columns[adata.obs.isnull().any()].tolist()
# print(f"Columns with NAs in the observation metadata (adata.obs): {na_columns_obs}")

# # Get a summary of missing values in each column
# na_summary = adata.obs.isnull().sum()
# na_summary = na_summary[na_summary > 0]
# print(f"Summary of NAs in each column:\n{na_summary}")

# pt_ids_with_na_pmi.unique()

# Define the mapping of Pt_ID to the new PMI values
# from: https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-023-00424-y/MediaObjects/43587_2023_424_MOESM1_ESM.pdf
pmi_replacement_values = {
    'D:1': 5.08,
    'D:10': 3.33,
    'D:3': 5.08,
    'D:15': 6.97
}

# Replace the NA values in the PMI column with the specified values
for pt_id, new_pmi in pmi_replacement_values.items():
    adata.obs.loc[adata.obs['Pt_ID'] == pt_id, 'PMI'] = new_pmi

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["Sex", "Race", "genotype_APOE", "CERAD"],
    continuous_covariate_keys=["percent.mito", "coded_Age", "nCount_RNA", "FreshBrainWeight", "PMI"],
    batch_key="SeqBatch"
)

model = scvi.model.SCVI(adata,
                        n_layers=2,
                        n_latent=30,
                        gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

model.save(dir_path="/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi2-model", 
           overwrite=True)

reserved_names = {'_index'}
if adata.raw is not None:
    raw_var_reserved = reserved_names.intersection(adata.raw.var.columns)
    print("Reserved names in raw.var:", raw_var_reserved)
    if '_index' in adata.raw.var.columns:
        adata.raw.var.rename(columns={'_index': 'index_raw_var'}, inplace=True)

adata.write("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi2_anndata.h5ad")

print("Done! :)")