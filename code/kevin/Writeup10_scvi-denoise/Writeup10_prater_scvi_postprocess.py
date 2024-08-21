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

adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scvi_anndata.h5ad")

# denoise the gene expression
# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
model = scvi.model.SCVI.load("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scvi-model", 
                             adata)

# adjust the donors covariates (leave the batch untouched)
# going based on https://github.com/linnykos/subject-de/blob/kevin/code/kevin/Writeup6_scvi/Writeup6_prater_scvi.py
# manually set the mode and mean
categorical_covariates = ["Sex", "Race", "APOEe4_status"]
# Store the original categories
original_categories = {covariate: adata.obs[covariate].cat.categories for covariate in categorical_covariates}

# Extracting the unique values of 'Pt_ID' from the obs attribute
unique_pt_ids = adata.obs["Pt_ID"].unique()
adata_obs_original = adata.obs.copy()

accumulator_matrix = np.zeros((adata.n_obs, adata.n_vars))
unique_batches = adata_obs_original["SeqBatch"].unique().tolist()

# Loop over each donor
for i in range(len(unique_pt_ids)):
    print(f"Working on subject {i + 1} out of {len(unique_pt_ids)}")
    
    pt = unique_pt_ids[i]
    subset = adata_obs_original[adata_obs_original['Pt_ID'] == pt]
    
    # grab the relevant values
    Sex_value = subset["Sex"].iloc[0]
    Race_value = subset["Race"].iloc[0]
    APOEe4_status_value = subset["APOEe4_status"].iloc[0]
    coded_Age_value = subset["coded_Age"].iloc[0]
    PMI_value = subset["PMI"].iloc[0]
    
    # set the values
    adata.obs["Sex"] = Sex_value
    adata.obs["Race"] = Race_value
    adata.obs["APOEe4_status"] = APOEe4_status_value
    adata.obs["coded_Age"] = coded_Age_value
    adata.obs["PMI"] = PMI_value
    
    for covariate in categorical_covariates:
        # Ensure the categorical variable type and categories are retained
        adata.obs[covariate] = adata.obs[covariate].astype('category')
        adata.obs[covariate] = adata.obs[covariate].cat.set_categories(original_categories[covariate])
    
    model.adata_manager.transfer_fields(adata_target=adata)
    denoised = model.get_normalized_expression(return_mean=False, 
                                               transform_batch=unique_batches,
                                               library_size=1000)
    denoised_array = denoised.to_numpy()
    
    accumulator_matrix += denoised_array

accumulator_matrix = accumulator_matrix / len(unique_pt_ids)
denoised_average = pd.DataFrame(accumulator_matrix)
denoised_average.index = adata.obs_names
denoised_average.columns = adata.var_names

# Assuming df is your DataFrame
denoised_average.to_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scvi.feather")
