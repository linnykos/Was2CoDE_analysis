# see Algorithm 1 in "Supervised principal component analysis: Visualization, classification andregression on subspaces and submanifolds"
# https://www.sciencedirect.com/science/article/pii/S0031320310005819

supervised_pca <- function(x, y, 
                           k = 2,
                           orthogonalize = TRUE,
                           scale_x = TRUE){
  stopifnot(is.matrix(x), is.matrix(y),
            nrow(x) == nrow(y))
  
  if(scale_x) x <- scale(x)
  n <- nrow(x)
  
  # H <- diag(n) - matrix(1/n, nrow = n, ncol = n)
  # half_mat <- crossprod(y, H %*% x)
  
  # see https://en.wikipedia.org/wiki/Centering_matrix
  x_centered <- scale(x, center = TRUE, scale = FALSE)
  half_mat <- crossprod(y, x_centered)
  
  Q <- crossprod(half_mat)
  eigen_res <- eigen(Q)
  U <- Re(eigen_res$vectors[,1:k])
  
  res <- x %*% U
  
  if(orthogonalize){
    svd_res <- svd(res)
    rotation_mat <- svd_res$v
    U <- U %*% rotation_mat
    res <- x %*% U
  }
  
  colnames(res) <- paste0("SPCA_", 1:k)
  colnames(U) <- paste0("SPCA_", 1:k)
  rownames(U) <- colnames(x)
  
  list(dimred = res, U = U)
}

form_onehot_classification_mat <- function(y){
  uniq_val <- sort(unique(y))
  k <- length(uniq_val)
  n <- length(y)
  
  mat <- matrix(0, nrow = n, ncol = k)
  for(j in 1:k){
    mat[which(y == uniq_val[j]),j] <- 1
  }
  colnames(mat) <- uniq_val
  if(length(names(y)) > 0) rownames(mat) <- names(y)
  
  mat
}

# based on https://cseweb.ucsd.edu/~dasgupta/251u/dist-embeddings-handout.pdf
distance_to_kernel <- function(mat){
  stopifnot(nrow(mat) == ncol(mat))
  n <- nrow(mat)
  
  H <- diag(n) - tcrossprod(rep(1/n, n))
  - (1/2) * H %*% mat %*% H
}