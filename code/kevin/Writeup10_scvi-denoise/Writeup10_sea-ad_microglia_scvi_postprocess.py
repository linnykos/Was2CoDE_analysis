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

adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scvi_anndata.h5ad")

# denoise the gene expression
# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
model = scvi.model.SCVI.load("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scvi-model", 
                             adata)

# adjust the donors covariates (leave the batch untouched)
# going based on https://github.com/linnykos/subject-de/blob/kevin/code/kevin/Writeup6_scvi/Writeup6_prater_scvi.py
# manually set the mode and mean
categorical_covariates = ["sex", "self_reported_ethnicity", "APOE4status"]
# Store the original categories
original_categories = {covariate: adata.obs[covariate].cat.categories for covariate in categorical_covariates}

unique_pt_ids = adata.obs['Pt_ID'].unique()

# Create an empty matrix of all 0's

# Loop over each unique donor

# Find

# Calculate mode for categorical covariates
for covariate in categorical_covariates:
    mode_value = adata.obs[covariate].mode()[0]
    adata.obs[covariate] = mode_value
    # Ensure the categorical variable type and categories are retained
    adata.obs[covariate] = adata.obs[covariate].astype('category')
    adata.obs[covariate] = adata.obs[covariate].cat.set_categories(original_categories[covariate])

# Calculate mean for continuous covariates
continuous_covariates = ["Ageatdeath", "PMI"]
for covariate in continuous_covariates:
    mean_value = adata.obs[covariate].mean()
    adata.obs[covariate] = mean_value

# ## double-checking
# adata.obs["Sex"].value_counts()
# adata.obs["coded_Age"].value_counts()

# push it through the decoder

# update the model's adata object

model.adata_manager.transfer_fields(adata_target=adata)
model.adata_manager.view_registry()

# push it through the decoder

unique_batches = adata.obs["SeqBatch"].unique().tolist()
denoised = model.get_normalized_expression(return_mean=False, 
                                           transform_batch=unique_batches,
                                           library_size=1000)

# Assuming df is your DataFrame
denoised.to_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scvi.feather")
