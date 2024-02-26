rm(list=ls())
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching.RData")
source("matching.R")

specific_genes <- c("ZDHHC21", "ADGRD1", "BPGM", "CYCS")

ntrials <- 1000
neuropath <- df[,which(colnames(df) != "cognition")]
lambda <- sum(eigen(stats::cov(neuropath))$values)
n <- nrow(neuropath)
gamma_seq <- seq(0, 1000, length.out = 5)

res_df <- res_df[specific_genes,]
sign_vec <- res_df[,"side"]

for(gamma in gamma_seq){
  pvalue_mat <- matrix(NA, ncol = length(specific_genes), nrow = ntrials)
  colnames(pvalue_mat) <- specific_genes
  
  for(trial in 1:ntrials){
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
                       gene = gene,
                       avg_mat = avg_mat,
                       signal_matches = matching_res$signal_matches)
    })
    
    res_df <- data.frame(
      pval = as.numeric(res_mat["pval",]),
      side = res_mat["side",]
    )
    rownames(res_df) <- colnames(res_mat)
    
    tmp <- sapply(1:nrow(res_df), function(j){
      if(res_df[j,"side"] != sign_vec[j]) return(1)
      res_df[j,"pval"]
    })
    names(tmp) <- specific_genes
    
    pvalue_mat[trial,names(tmp)] <- tmp
  }
  
  pvalue_vec <- apply(pvalue_mat, 2, max)
}
