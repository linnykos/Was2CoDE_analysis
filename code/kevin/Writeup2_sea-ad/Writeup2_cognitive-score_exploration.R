rm(list=ls())
library(openxlsx)

metadata <- openxlsx::read.xlsx("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_cohort_donor_metadata_082222.xlsx")
rownames(metadata) <- metadata$Donor.ID
harmonized_scores <- openxlsx::read.xlsx("~/Dropbox/Collaboration-and-People/alzheimer/data/joey-harmonized-cognitive/sea_ad_cognitive_slopes_79.xlsx")
rownames(harmonized_scores) <- harmonized_scores$donor_name
neuropath <- read.csv("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv")
rownames(neuropath) <- neuropath$Donor.ID

path_idx <- grep("^percent.*area_Grey", colnames(neuropath))
colnames(neuropath)[path_idx]
neuropath <- neuropath[,path_idx]
summary(harmonized_scores)

included_idx <- which(metadata$Donor.ID %in% harmonized_scores$donor_name)
metadata <- metadata[included_idx,]

ad_idx <- intersect(
  which(metadata$Braak %in% c("Braak V", "Braak VI")),
  which(metadata$CERAD.score %in% c("Moderate", "Frequent"))
)
metadata <- metadata[ad_idx,]
harmonized_scores <- harmonized_scores[rownames(metadata),]
neuropath <- neuropath[rownames(metadata),]

summary(neuropath)
summary(harmonized_scores)
hist(harmonized_scores$slope_zmem0)

###################

neuropath2 <- scale(neuropath)
pca_res <- stats::prcomp(neuropath2)

col_palette <- viridis::viridis(10)
val_range <- range(harmonized_scores$slope_zmem0)
val_breaks <- seq(val_range[1], val_range[2], length.out = 10)
col_vec <- sapply(harmonized_scores$slope_zmem0, function(x){
  col_palette[which.min(abs(x - val_breaks))]
})

med_val <- stats::median(harmonized_scores$slope_zmem0)
col_vec2 <- sapply(harmonized_scores$slope_zmem0, function(x){
  if(x <= med_val) 2 else 4
})

par(mfrow = c(1,2), mar = c(4,4,4,0.5))
plot(pca_res$x[,1], pca_res$x[,2],
     pch = 16, 
     asp = T, 
     cex = 2,
     xlab = "PCA1 from neuropath",
     ylab = "PCA2 from neuropath")
points(pca_res$x[,1], pca_res$x[,2],
       pch = 16, 
       col = col_vec,
       cex = 1.5)
plot(pca_res$x[,1], pca_res$x[,2],
     pch = 16, 
     asp = T, 
     cex = 2,
     xlab = "PCA1 from neuropath",
     ylab = "PCA2 from neuropath")
points(pca_res$x[,1], pca_res$x[,2],
       pch = 16, 
       col = col_vec2,
       cex = 1.5)

colname_vec <- sapply(colnames(neuropath), function(x){
  strsplit(x, split = "\\.")[[1]][2]
})
col_palette3 <- grDevices::colorRampPalette(c("purple", "gray", "orange"))(10)
par(mfrow = c(2,4), mar = c(0.5, 0.5, 3, 0.5))
for(j in 1:7){
  val_range <- range(neuropath2[,j])
  val_breaks <- seq(val_range[1], val_range[2], length.out = 10)
  col_vec3 <- sapply(neuropath2[,j], function(x){
    col_palette3[which.min(abs(x - val_breaks))]
  })
  
  plot(pca_res$x[,1], pca_res$x[,2],
       pch = 16, 
       asp = T, 
       cex = 2,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       main = colname_vec[j])
  points(pca_res$x[,1], pca_res$x[,2],
         pch = 16, 
         col = col_vec3,
         cex = 1.5)
}



###############

set.seed(10)
isomap_res <- dimRed::embed(neuropath2, "Isomap", knn = 5)
isomap_mat <- isomap_res@data@data

par(mfrow = c(1,2), mar = c(4,4,4,0.5))
plot(isomap_mat[,1], isomap_mat[,2],
     pch = 16, 
     asp = T, 
     cex = 2,
     xlab = "Isomap1 from neuropath",
     ylab = "Isomap2 from neuropath")
points(isomap_mat[,1], isomap_mat[,2],
       pch = 16, 
       col = col_vec,
       cex = 1.5)
plot(isomap_mat[,1], isomap_mat[,2],
     pch = 16, 
     asp = T, 
     cex = 2,
     xlab = "Isomap1 from neuropath",
     ylab = "Isomap2 from neuropath")
points(isomap_mat[,1], isomap_mat[,2],
       pch = 16, 
       col = col_vec2,
       cex = 1.5)

colname_vec <- sapply(colnames(neuropath), function(x){
  strsplit(x, split = "\\.")[[1]][2]
})
col_palette3 <- grDevices::colorRampPalette(c("purple", "gray", "orange"))(10)
par(mfrow = c(2,4), mar = c(0.5, 0.5, 3, 0.5))
for(j in 1:7){
  val_range <- range(neuropath2[,j])
  val_breaks <- seq(val_range[1], val_range[2], length.out = 10)
  col_vec3 <- sapply(neuropath2[,j], function(x){
    col_palette3[which.min(abs(x - val_breaks))]
  })
  
  plot(isomap_mat[,1], isomap_mat[,2],
       pch = 16, 
       asp = T, 
       cex = 2,
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       main = colname_vec[j])
  points(isomap_mat[,1], isomap_mat[,2],
         pch = 16, 
         col = col_vec3,
         cex = 1.5)
}

