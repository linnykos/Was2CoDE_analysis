rm(list=ls())
source("kevin/Writeup7_decomp-test/compositional_test.R")

m <- 20 # number of donors in each class
n_each <- 200 # number of cells per donor

set.seed(10)
x1 <- lapply(1:m, function(ii){
  stats::rnorm(n_each)
})
# x2 <- lapply(1:m, function(ii){
#   c(stats::rnorm(round(n_each/2), mean = -2),
#     stats::rnorm(round(n_each/2), mean = 2))
# })
# x2 <- lapply(1:m, function(ii){
#   stats::rnorm(n_each, sd = 1)
# })
# x2 <- lapply(1:m, function(ii){
#   stats::rnorm(n_each, sd = 10)
# })
x2 <- lapply(1:m, function(ii){
  stats::rnorm(n_each, sd = 1)
})
x <- c(unlist(x1), unlist(x2))
df <- data.frame(x)
colnames(df) <- "value"
df$donor <- rep(paste0("pt", 1:(2*m)), each = n_each)
df$cc <- rep(paste0("g", 1:2), each = m*n_each)

compositional_test(df,
                   value_col = "value",
                   donor_col = "donor",
                   cc_col = "cc",
                   minPts = sqrt(nrow(df)),
                   max_cluster = 6)
