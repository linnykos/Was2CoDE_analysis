matching_func <- function(df,
                          class_covariate){
  stopifnot(class_covariate %in% colnames(df))
  idx <- which(colnames(df) == class_covariate)
  colnames(df)[idx] <- "zz"
  
  set.seed(10)
  fm_res <- optmatch::fullmatch(zz ~ ., data = df)
  # print(fm_res, grouped = TRUE) # what a garbage method... I can't store this?
  signal_matches <- sapply(unique(fm_res), function(x){
    names(fm_res)[which(fm_res == x)[1:2]] # TODO: we're throw away extra donors for now
  })
  # make sure the order in each column is correct. cognition=TRUE on top
  for(j in 1:ncol(signal_matches)){
    vec <- signal_matches[,j]
    if(df[vec[1],"zz"] == FALSE) vec <- vec[c(2,1)]
    signal_matches[,j] <- vec
  }
  
  neuropath <- df[,which(colnames(df) != "zz")]
  neuropath_false <- neuropath[!df$zz,]
  nn_res_false <- RANN::nn2(neuropath_false, k = 2)$nn.idx
  rownames(nn_res_false) <- rownames(neuropath_false)
  nn_res_false <- rownames(nn_res_false)[nn_res_false[,2]]
  names(nn_res_false) <- rownames(neuropath_false) # let's just ignore doubling for now
  
  neuropath_true <- neuropath[df$zz,]
  nn_res_true <- RANN::nn2(neuropath_true, k = 2)$nn.idx
  rownames(nn_res_true) <- rownames(neuropath_true)
  nn_res_true <- rownames(nn_res_true)[nn_res_true[,2]]
  names(nn_res_true) <- rownames(neuropath_true) # let's just ignore doubling for now
  
  bg_matches <- c(nn_res_false, nn_res_true)
  
  list(signal_matches = signal_matches,
       bg_matches = bg_matches)
}

compute_wilcoxon <- function(bg_matches,
                             gene,
                             avg_mat,
                             signal_matches){
  # compute mean expression for each donor
  expression_vec <- avg_mat[,gene]
  names(expression_vec) <- rownames(avg_mat)
  
  # compute the signal pairings
  signal_vec <- sapply(1:ncol(signal_matches), function(j){
    donor1 <- signal_matches[1,j]
    donor2 <- signal_matches[2,j]
    
    expression_vec[donor1] - expression_vec[donor2]
  })
  
  # compute the bg pairings
  bg_vec <- sapply(1:length(bg_matches), function(j){
    donor1 <- names(bg_matches)[j]
    donor2 <- bg_matches[j]
    
    expression_vec[donor1] - expression_vec[donor2]
  })
  bg_vec <- c(bg_vec, -bg_vec)
  bg_vec <- bg_vec[!duplicated(bg_vec)]
  
  wilcox_twosided <- stats::wilcox.test(
    x = signal_vec,
    y = bg_vec,
    alternative = "two.sided"
  )$p.value
  wilcox_less <- stats::wilcox.test(
    x = signal_vec,
    y = bg_vec,
    alternative = "less"
  )$p.value
  wilcox_greater <- stats::wilcox.test(
    x = signal_vec,
    y = bg_vec,
    alternative = "greater"
  )$p.value
  
  c(twosided = wilcox_twosided,
    less = wilcox_less,
    greater = wilcox_greater)
  
}