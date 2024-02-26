rm(list=ls())
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching.RData")
source("matching.R")

# specific_genes <- c("ZDHHC21", "ADGRD1", "BPGM", "CYCS")

gene_names <- rownames(res_df)
sun_sheet <- openxlsx::read.xlsx(
  xlsxFile = "~/kzlinlab/projects/subject-de/data/1-s2.0-S0092867423009716-mmc1.xlsx",
  sheet = "Page 10.DEGs_AD"
)
sun_genes <- sort(unique(sun_sheet[which(sun_sheet[,"fdr"] <= 0.05),"row.names"]))
sun_genes <- intersect(gene_names, sun_genes)

hk_df <- read.csv("~/kzlinlab/projects/subject-de/data/Housekeeping_GenesHuman.csv", sep = ";")
hk_genes <- hk_df[,"Gene.name"]
hk_genes <- intersect(gene_names, hk_genes)

sun_genes <- sun_genes[order(abs(abs(res_df[sun_genes,"test_stat"]) - 0.5), decreasing = F)[1:50]]
hk_genes <- hk_genes[order(abs(abs(res_df[hk_genes,"test_stat"]) - 0.5), decreasing = F)[1:50]]

ntrials <- 5000
neuropath <- df[,which(colnames(df) != "cognition")]
n <- nrow(neuropath)
gamma_seq <- round(10^(seq(-2, 1, length.out = 10)),2)

specific_genes <- unique(c(sun_genes, hk_genes))
res_df_target <- res_df[specific_genes,]
class_covariate <- "cognition"
other_var_idx <- which(colnames(df) != class_covariate)

sa_mat <- matrix(NA, nrow = length(specific_genes), ncol = length(gamma_seq))
rownames(sa_mat) <- specific_genes
colnames(sa_mat) <- paste0("gamma:", gamma_seq)

for(kk in 1:length(gamma_seq)){
  gamma <- gamma_seq[kk]
  print(paste0("Working on gamma: ", gamma))
  teststat_mat <- matrix(NA, ncol = length(specific_genes), nrow = ntrials)
  colnames(teststat_mat) <- specific_genes
  
  for(trial in 1:ntrials){
    if(trial %% floor(ntrials/10) == 0) cat('*')
    
    df2 <- df
    set.seed(trial)
    for(j in other_var_idx){
      sd_val <- stats::sd(df[,j])
      df2[,j] <- df2[,j] + stats::rnorm(n, mean = 0, sd = gamma*sd_val)
      # df2[,j] <- scale(df2[,j])
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
  names(teststat_vec) <-  colnames(teststat_mat)
  
  sa_mat[ names(teststat_vec),kk] <- teststat_vec
  
  # round(cbind(res_df_target[,"test_stat"], teststat_vec),2)
}

sa_mat <- cbind(res_df_target[,"test_stat"], sa_mat)
# sa_mat_safe <- sa_mat
# 
# round(sa_mat,2)
# for(i in 1:nrow(sa_mat)){
#   sa_mat[i,] <- sort(sa_mat[i,],decreasing = F)
# }
# 
# steady_genes <- c("TCFL5", "XRCC2", "SRSF3", "TCP1", "SLC41A2", "EFCAB6")
# steady_genes[steady_genes %in% sun_genes]
# # XRCC2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8079392/
# # EFCAB6: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10532410/
# 
# fast_genes <- c("SOX8", "BCL11B", "TNS1", "POLR2I", "PTGES2", "BAP1", "PRMT1")
# fast_genes[steady_genes %in% hk_genes]
# # TNS1 and BAP1

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(sa_mat_safe,
     date_of_run, session_info, 
     file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching_sensitivity.RData")
