import os
import tempfile
import scanpy as sc
import scvi
import seaborn as sns
import torch
import anndata
import pandas as pd

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = "/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/"

adata = anndata.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia.h5ad")
adata

sc.pp.filter_genes(adata, min_counts=3)
adata.layers["counts"] = adata.X.copy()  # preserve counts

# layer_data = adata.layers['counts']
# sub_matrix = pd.DataFrame(layer_data[:20, :20].toarray() if hasattr(layer_data, "toarray") else layer_data[:20, :20],
#                           index=adata.obs_names[:20],
#                           columns=adata.var_names[:20])
# sub_matrix

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)
adata

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["donor_id", "sex", "self_reported_ethnicity", "Lewy_body_disease_pathology", "APOE4_status"],
    continuous_covariate_keys=["AgeAtDeath", "YearsOfEducation", "FractionMitochrondrialUMIs", "percent_AT8", "percent_pTDP43", 
                               "percent_Iba1", "percent_6e10", "percent_aSyn", "percent_NeuN", "percent_GFAP"]
)

model = scvi.model.SCVI(adata,
                        n_latent=15)
model

model.train()

model_dir = os.path.join("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/scvi_model")
model.save(model_dir, overwrite=True)

print("Done! :)")