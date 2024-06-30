import scanpy as sc
import scipy 

adata = sc.read_h5ad("/home/users/kzlin/kzlinlab/data/sea-ad/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad")
print(adata)
print("cell_prep_type")
print(adata.obs['cell_prep_type'].value_counts())

print("sample_name")
print(adata.obs['sample_name'].value_counts())

print("method")
print(adata.obs['method'].value_counts())

print("library_prep")
print(adata.obs['library_prep'].value_counts())

print("Supertype")
print(adata.obs['Supertype'].value_counts())

print("Subclass")
print(adata.obs['Subclass'].value_counts())

print("Cell names")
print(adata.obs_names[0:10])

print("Feature names")
print(adata.var_names[0:10])

# Create a boolean mask for cells with '10xMulti' in 'method' and 'Microglia-PVM' in 'Subclass'
mask = (adata.obs['method'] == '10xMulti') & (adata.obs['Subclass'] == 'Microglia-PVM')

# Subset the AnnData object
adata_subset = adata[mask]

print("Cell subset names")
print(adata_subset.obs_names[0:10])

# Save the subsetted AnnData object to an HDF5 file
adata_subset.write_h5ad('/home/users/kzlin/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_SEAAD_MTG_ATACseq_10xMulti_Microglia-PVM.h5ad')