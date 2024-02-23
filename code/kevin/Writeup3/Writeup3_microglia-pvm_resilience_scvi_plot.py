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
sc.settings.figdir = "/home/users/kzlin/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/"

adata = anndata.read_h5ad("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia.h5ad")
sc.pp.filter_genes(adata, min_counts=3)
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)
adata

model = scvi.model.SCVI.load("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/scvi_model", 
                            adata=adata)
model

SCVI_LATENT_KEY = "X_scVI"
latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

print("Saving peakVI as CSV")
cell_names = adata.obs_names
df = pd.DataFrame(latent, index = cell_names)
df.to_csv("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_microglia-pvm_resilience_scvi.csv")

SCVI_NORMALIZED_KEY = "scvi_normalized"
gene_names = adata.var_names
adata.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)
normalized_mat = adata.layers[SCVI_NORMALIZED_KEY]
df = pd.DataFrame(normalized_mat, index = cell_names, columns = gene_names)
df.to_csv("/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_microglia-pvm_resilience_scvi_normalized-mat.csv")

# run PCA then generate UMAP plots
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)

filename = "scvi_donor-id.png"
sc.pl.umap(
    adata,
    color=["donor_id"],
    frameon=False,
    save=filename
)

filename = "scvi_resilient_sex_age.png"
sc.pl.umap(
    adata,
    color=["resilient", "sex", "AgeAtDeath"],
    ncols=3,
    frameon=False,
    save=filename
)