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

adata = ad.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup11/Writeup11_rosmap_scvi_anndata.h5ad")

# denoise the gene expression
# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html
model = scvi.model.SCVI.load("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup11/Writeup11_rosmap_scvi-model", 
                             adata)

# adjust the donors covariates (leave the batch untouched)
# going based on https://github.com/linnykos/subject-de/blob/kevin/code/kevin/Writeup6_scvi/Writeup6_prater_scvi.py
# manually set the mode and mean
categorical_covariates = ["sex", "race", "APOEe4_status"]
# Store the original categories
original_categories = {covariate: adata.obs[covariate].cat.categories for covariate in categorical_covariates}

# Extracting the unique values of 'donor_id' from the obs attribute
unique_pt_ids = adata.obs["donor_id"].unique()
adata_obs_original = adata.obs.copy()

accumulator_matrix = np.zeros((adata.n_obs, adata.n_vars))
unique_batches = adata_obs_original["assay"].unique().tolist()

# Loop over each donor
for i in range(len(unique_pt_ids)):
    print(f"Working on subject {i + 1} out of {len(unique_pt_ids)}")
    
    pt = unique_pt_ids[i]
    subset = adata_obs_original[adata_obs_original['donor_id'] == pt]
    
    # grab the relevant values
    sex_value = subset["sex"].iloc[0]
    ethnicity_value = subset["race"].iloc[0]
    APOE4status_value = subset["APOEe4_status"].iloc[0]
    Ageatdeath_value = subset["age_death"].iloc[0]
    PMI_value = subset["pmi"].iloc[0]
    
    # set the values
    adata.obs["sex"] = sex_value
    adata.obs["race"] = ethnicity_value
    adata.obs["APOEe4_status"] = APOE4status_value
    adata.obs["age_death"] = Ageatdeath_value
    adata.obs["pmi"] = PMI_value
    
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
denoised_average.to_feather("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup11/Writeup11_rosmap_scvi.feather")
