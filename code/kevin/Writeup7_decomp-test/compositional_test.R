library(tidyr)
library(dplyr)
library(compositions)
library(dbscan)
library(RNOmni)

compositional_test <- function(df,
                               value_col = "value",
                               donor_col = "donor",
                               cc_col = "cc",
                               minPts = sqrt(nrow(df)),
                               max_cluster = 6){

  colnames(df)[which(colnames(df) == donor_col)] <- "donor"
  colnames(df)[which(colnames(df) == cc_col)] <- "cc"
  
  tmp_mat <- matrix(df[,value_col], ncol = 1)
  hclust_res <- dbscan::hdbscan(tmp_mat,
                                minPts = minPts)
  cluster_res <- hclust_res$cluster
  
  # handle no clustering
  if(length(unique(cluster_res)) == 1 & all(cluster_res == 0)){
    hclust_res <- dbscan::hdbscan(matrix(df[,value_col], ncol = 1), 
                                  minPts = 5)
  }
  
  if(length(unique(cluster_res)) == 1) return(1)
  
  # handle too many clusters
  if(length(unique(cluster_res)) > max_cluster+1){
    print("Too many clusters...")
    print(table(cluster_res))
    
    ## TODO: This needs to be tested...
    cluster_res <- stats::cutree(hclust_res$hc, 
                                 k = max_cluster)
  }
  
  # assign noise points to clusters
  n <- length(cluster_res)
  nonzero_idx <- which(cluster_res != 0)
  nonzero_cluster <- cluster_res[nonzero_idx]
  if(length(nonzero_idx) != n){
    for(i in setdiff(1:n, nonzero_idx)){
      val <- df[i,value_col]
      cluster_res[i] <- nonzero_cluster[which.min(abs(df[nonzero_idx,value_col] - val))]
    }
  }
  
  cluster_res <- paste0("c", cluster_res)
  print("Final clustering")
  print(table(cluster_res))
  
  # convert to data frame
  df$cluster <- cluster_res
  df2 <- tidyr::as_tibble(df[,c("cluster", "donor", "cc")]) %>% 
    dplyr::group_by(donor, cluster, cc) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = cluster, values_from = count, values_fill = 0)
  
  df2 <- as.data.frame(df2)
  
  # perform compositional analysis
  xc <- compositions::acomp(df2[, sort(unique(df$cluster))])
  xc <- compositions::zeroreplace(xc)
  
  df2$comp <- xc
  codalm <- stats::lm(compositions::clr(comp)~cc, data=df2) 
  tmp <- summary(codalm)
  pvalue_vec <- sapply(tmp, function(obj){
    obj$coefficients[2,"Pr(>|t|)"]
  })
  
  RNOmni::OmniP(pvalue_vec)
}