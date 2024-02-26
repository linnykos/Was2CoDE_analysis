rm(list=ls())
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching.RData")
source("matching.R")

specific_genes <- c("ZDHHC21", "ADGRD1", "BPGM", "CYCS")

ntrials <- 1000
neuropath <- df[,which(colnames(df) != "cognition")]
n <- nrow(neuropath)
gamma_seq <- seq(0, 1000, length.out = 5)

res_df_target <- res_df[specific_genes,]
class_covariate <- "cognition"
other_var_idx <- which(colnames(df) != class_covariate)

for(gamma in gamma_seq){
  teststat_mat <- matrix(NA, ncol = length(specific_genes), nrow = ntrials)
  colnames(teststat_mat) <- specific_genes
  
  for(trial in 1:ntrials){
    if(trial %% floor(ntrials/10) == 0) cat('*')
    
    df2 <- df
    set.seed(trial)
    for(j in other_var_idx){
      sd_val <- stats::sd(df[,j])
      df2[,j] <- df2[,j] + stats::rnorm(n, mean = 0, sd = gamma*sd_val)
    }
    
    matching_res <- matching_func(
      df = df2,
      class_covariate = "cognition"
    )
    
    res_mat <- sapply(specific_genes, function(gene){
      compute_test(bg_matches = matching_res$bg_matches,
                       gene = gene,
                       avg_mat = avg_mat,
                       signal_matches = matching_res$signal_matches)
    })
    
    res_df <- data.frame(t(res_mat))
    rownames(res_df) <- colnames(res_mat)
    
    teststat_mat[trial,rownames(res_df)] <- res_df$test_stat
  }
  
  teststat_vec <- sapply(1:ncol(teststat_mat), function(j){
    target_val <- res_df_target[j,"test_stat"]
    tmp <- teststat_mat[,j]
    
    if(target_val > 0){
      tmp <- pmin(tmp, target_val)
    } else {
      tmp <- pmax(tmp, target_val)
    }
    
    tmp[which.max(abs(tmp - target_val))]
  })
}
