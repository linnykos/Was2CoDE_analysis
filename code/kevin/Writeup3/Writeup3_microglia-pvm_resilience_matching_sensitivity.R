rm(list=ls())
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching.RData")
source("matching.R")

ntrials <- 1000
neuropath <- df[,which(colnames(df) != "cognition")]
lambda <- eigen(stats::cov(neuropath))$values[1]
n <- nrow(neuropath)
gamma_seq <- seq(0, 1000, length.out = 5)

res_mat2 <- res_mat[,specific_genes]
sign_vec <- sapply(specific_genes, function(gene){
  if(res_mat2["less",gene] < 0.01) "less" else "greater"
})

for(gamma in gamma_seq){
  pvalue_mat <- matrix(NA, ncol = length(specific_genes), nrow = ntrials)
  colnames(pvalue_mat) <- specific_genes
  
  for(trial in 1:ntrials){
    print(trial)
    if(trial %% floor(ntrials/10) == 0) cat('*')
    
    set.seed(trial)
    unobs <- stats::rnorm(n, mean = 0)
    unobs <- scale(unobs)*sqrt(gamma*lambda)
    
    df2 <- df; df2$unobs <- unobs
    
    matching_res <- matching_func(
      df = df2,
      class_covariate = "cognition"
    )
    
    res_mat <- sapply(specific_genes, function(gene){
      compute_wilcoxon(bg_matches = matching_res$bg_matches,
                       donor_list = donor_list,
                       gene = gene,
                       mat = mat,
                       signal_matches = matching_res$signal_matches)
    })
    
    tmp <- sapply(1:ncol(res_mat), function(j){
      sign_true <- sign_vec[j]
      sign_opposite <- setdiff(c("less", "greater"), sign_true)
      
      if(res_mat[sign_true,j] > res_mat[sign_opposite,j]) return(1)
      res_mat["twosided",j]
    })
    names(tmp) <- specific_genes
    
    pvalue_mat[trial,names(tmp)] <- tmp
  }
  
  pvalue_vec <- apply(pvalue_mat, 2, max)
}
