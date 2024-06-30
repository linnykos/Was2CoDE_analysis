rm(list=ls())

library(meanShiftR)
library(tidyr)
library(dplyr)
library(compositions)

set.seed(10)
x <- matrix(runif(20),10,2)
classification <- meanShiftR::meanShift(queryData = x,
                                        trainData = x)
classification


set.seed(10)
x <- matrix(c(
  stats::rnorm(100),
  stats::rnorm(100, mean = 2)
), ncol = 1)

# bw <- stats::bw.nrd(x) # https://stat.ethz.ch/R-manual/R-devel/library/stats/html/bandwidth.html
# # bw <- KernSmooth::dpik(x)
# classification <- meanShiftR::meanShift(queryData = x,
#                                         trainData = x,
#                                         bandwidth = bw)
# table(classification$assignment)

set.seed(10)
hclust_res <- dbscan::hdbscan(x, minPts = sqrt(length(x)))
cluster_res <- hclust_res$cluster
# I'm not sure how the "predict" function works

# assign cluster 0's to the closest cluster
n <- length(cluster_res)
nonzero_idx <- which(cluster_res != 0)
nonzero_cluster <- cluster_res[nonzero_idx]
if(length(nonzero_idx) != n){
  for(i in setdiff(1:n, nonzero_idx)){
    val <- x[i,1]
    cluster_res[i] <- nonzero_cluster[which.min(abs(x[nonzero_idx,1] - val))]
  }
}
table(cluster_res)

par(mfrow = c(1,2))
plot(x[,1], pch = 16, col = cluster_res)
hist(x); for(i in 1:n){rug(x[i,1], col = cluster_res[i])}

#################

# https://github.com/kharchenkolab/cacoa/blob/main/R/coda.R?

m <- 20 # number of donors in each class
n_each <- 200 # number of cells per donor

x1 <- lapply(1:m, function(ii){
  stats::rnorm(n_each)
})
x2 <- lapply(1:m, function(ii){
  c(stats::rnorm(round(n_each/2), mean = -2),
    stats::rnorm(round(n_each/2), mean = 2))
})
x <- c(unlist(x1), unlist(x2))
x <- matrix(x, ncol = 1)

set.seed(10)
hclust_res <- dbscan::hdbscan(x, minPts = round(sqrt(length(x))))
cluster_res <- hclust_res$cluster
table(cluster_res)

# assign cluster 0's to the closest cluster
n <- length(cluster_res)
nonzero_idx <- which(cluster_res != 0)
nonzero_cluster <- cluster_res[nonzero_idx]
if(length(nonzero_idx) != n){
  for(i in setdiff(1:n, nonzero_idx)){
    val <- x[i,1]
    cluster_res[i] <- nonzero_cluster[which.min(abs(x[nonzero_idx,1] - val))]
  }
}
table(cluster_res)

par(mfrow = c(1,2))
plot(x[,1], pch = 16, col = cluster_res)
hist(x); for(i in 1:n){rug(x[i,1], col = cluster_res[i])}

df <- data.frame(cbind(x, cluster_res))
colnames(df) <- c("value", "cluster")
df$cluster <- paste0("c", df$cluster)
df$donor <- rep(paste0("pt", 1:(2*m)), each = n_each)
df$cc <- rep(paste0("g", 1:2), each = m*n_each)

df2 <- tidyr::as_tibble(df[,c("cluster", "donor", "cc")]) %>% 
  dplyr::group_by(donor, cluster, cc) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = cluster, values_from = count, values_fill = 0)

df2 <- as.data.frame(df2)

xc <- compositions::acomp(df2[, c("c1", "c2", "c3", "c4", "c5", "c6")])
xc <- compositions::zeroreplace(xc)
# See Section 7.2.1.3 in http://ndl.ethernet.edu.et/bitstream/123456789/39269/1/K.%20Gerald%20van%20den.pdf
summary(xc)

# https://search.r-project.org/CRAN/refmans/compositions/html/compmiss.html
# to handle "BDL": Below detection limit


df2$comp <- xc
codalm = lm(compositions::alr(comp)~cc, data=df2) 
tmp <- summary(codalm)
pvalue_vec <- sapply(tmp, function(obj){
  obj$coefficients[2,"Pr(>|t|)"]
})

# maybe output the omnibus p-value?
library(RNOmni)
RNOmni::OmniP(pvalue_vec)
