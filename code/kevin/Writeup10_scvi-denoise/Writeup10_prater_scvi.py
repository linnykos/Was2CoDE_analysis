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

print("Last run with scvi-tools version:", scvi.__version__)

file_path = "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_cleaned.h5ad"
adata = ad.read_h5ad(file_path)

adata.obs["SeqBatch"] = adata.obs["SeqBatch"].astype('category')
adata.obs["Sex"] = adata.obs["Sex"].astype('category')
adata.obs["Race"] = adata.obs["Race"].astype('category')
adata.obs["APOEe4_status"] = adata.obs["APOEe4_status"].astype('category')

# adata.obs["Race"].value_counts()

# adding raw counts for referring to it in the future
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000,
    batch_key="Pt_ID",
    subset=True
)

adata

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["Sex", "Race", "APOEe4_status"],
    continuous_covariate_keys=["coded_Age", "PMI"],
    batch_key="SeqBatch"
)

model = scvi.model.SCVI(adata,
                        n_layers=2,
                        n_latent=30,
                        gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

model.save(dir_path="/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scvi-model", 
           overwrite=True)

reserved_names = {'_index'}
if adata.raw is not None:
    raw_var_reserved = reserved_names.intersection(adata.raw.var.columns)
    print("Reserved names in raw.var:", raw_var_reserved)
    if '_index' in adata.raw.var.columns:
        adata.raw.var.rename(columns={'_index': 'index_raw_var'}, inplace=True)

adata.write("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scvi_anndata.h5ad")

print("Done! :)")