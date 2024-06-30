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

count_matrix = adata.X
if scipy.sparse.issparse(count_matrix):
    count_matrix = count_matrix.todense()

# Print the first 10-by-10 submatrix
submatrix = count_matrix[:10, :10]
print("Submatrix")
print(submatrix)