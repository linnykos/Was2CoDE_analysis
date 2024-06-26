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

adata.obs["Pt_ID"] = adata.obs["Pt_ID"].astype('category')

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

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="Pt_ID"
)

model = scvi.model.SCVI(adata,
                        n_layers=2,
                        n_latent=30,
                        gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

model.save(dir_path="/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-noCovariates-model", 
           overwrite=True)

reserved_names = {'_index'}
if adata.raw is not None:
    raw_var_reserved = reserved_names.intersection(adata.raw.var.columns)
    print("Reserved names in raw.var:", raw_var_reserved)
    if '_index' in adata.raw.var.columns:
        adata.raw.var.rename(columns={'_index': 'index_raw_var'}, inplace=True)

adata.write("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-noCovariates_anndata.h5ad")

print("Done! :)")