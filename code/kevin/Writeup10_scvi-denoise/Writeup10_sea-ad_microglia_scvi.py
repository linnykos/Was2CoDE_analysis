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

file_path = "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_cleaned.h5ad"
adata = ad.read_h5ad(file_path)

adata.obs["assay"] = adata.obs["assay"].astype('category')
adata.obs["sex"] = adata.obs["sex"].astype('category')
adata.obs["self_reported_ethnicity"] = adata.obs["self_reported_ethnicity"].astype('category')
adata.obs["APOE4status"] = adata.obs["APOE4status"].astype('category')
adata.obs["donor_id"] = adata.obs["donor_id"].astype('category')

# adata.obs["self_reported_ethnicity"].value_counts()
# adata.obs["donor_id"].value_counts()

# adding raw counts for referring to it in the future
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000,
    batch_key="donor_id",
    subset=True
)

adata

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["sex", "self_reported_ethnicity", "APOE4status"],
    continuous_covariate_keys=["Ageatdeath", "PMI"],
    batch_key="assay"
)

model = scvi.model.SCVI(adata,
                        n_layers=2,
                        n_latent=30,
                        gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

model.save(dir_path="/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scvi-model", 
           overwrite=True)

reserved_names = {'_index'}
if adata.raw is not None:
    raw_var_reserved = reserved_names.intersection(adata.raw.var.columns)
    print("Reserved names in raw.var:", raw_var_reserved)
    if '_index' in adata.raw.var.columns:
        adata.raw.var.rename(columns={'_index': 'index_raw_var'}, inplace=True)

adata.write("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scvi_anndata.h5ad")

print("Done! :)")