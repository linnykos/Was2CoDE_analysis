
gene_index <- 738  #gene index for "SMO"

# Extract gene data from the dist_list
gene_distance <- dist_list$distance[gene_index, ,]
gene_location <- dist_list$location[gene_index, ,]
gene_size <- dist_list$size[gene_index, ,]
gene_shape <- dist_list$shape[gene_index, ,] 

# Initialize matrices to store results
n <- nrow(meta_ind)  # Number of donors
location_ratio <- matrix(NA, n, n)
size_ratio <- matrix(NA, n, n)
shape_ratio <- matrix(NA, n, n)

# Calculate the ratios for each pair of donors
for (i in 1:n) {
  for (j in 1:n) {
    # Calculate squared distance to avoid division by zero and for use in ratio calculations
    distance_squared <- gene_distance[i, j]^2
      location_diff <- gene_location[i,j]
      size_diff <- gene_size[i,j]
      shape_diff <- gene_shape[i,j] 
      
      location_ratio[i, j] <- location_diff / distance_squared
      size_ratio[i, j] <- size_diff / distance_squared
      shape_ratio[i, j] <- shape_diff / distance_squared
  }
}
rownames(location_ratio) <- meta_ind$individual
rownames(size_ratio) <- meta_ind$individual
rownames(shape_ratio) <- meta_ind$individual
colnames(location_ratio) <- meta_ind$individual
colnames(size_ratio) <- meta_ind$individual
colnames(shape_ratio) <- meta_ind$individual

dementia_indices <- which(meta_ind$CognitiveStatus == "Dementia")
no_dementia_indices <- which(meta_ind$CognitiveStatus == "No_dementia")

location_ratio_dnd <- location_ratio[dementia_indices, no_dementia_indices]
size_ratio_dnd <- size_ratio[dementia_indices, no_dementia_indices]
shape_ratio_dnd <- size_ratio[dementia_indices, no_dementia_indices]

location_ratio_dd <- location_ratio[dementia_indices, dementia_indices]
size_ratio_dd <- size_ratio[dementia_indices, dementia_indices]
shape_ratio_dd <- size_ratio[dementia_indices, dementia_indices]

location_ratio_ndnd <- location_ratio[no_dementia_indices, no_dementia_indices]
size_ratio_ndnd <- size_ratio[no_dementia_indices, no_dementia_indices]
shape_ratio_ndnd <- size_ratio[no_dementia_indices, no_dementia_indices]