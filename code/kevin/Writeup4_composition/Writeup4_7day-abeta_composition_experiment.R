rm(list=ls())
library(Seurat)
library(radEmu)

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_7day-abeta_simplified.RData")

set.seed(10)

# https://github.com/statdivlab/radEmu_supplementary/blob/main/fig456/0-process.R
# formatting the data
Y <- table(ss_data_norm$Sample_ID, ss_data_norm$SCT_snn_res.0.5)
Y <- data.frame(Y)

write.csv(Y, 
          row.names = TRUE,
          file = "~/kzlinlab/projects/subject-de/out/kevin/Writeup4/Writeup4_7day-abeta_composition.csv")

group_list <- list(
  het_7d_abeta = "het_7d_abeta_br1_tr2",
  het_ut = c("het_ut_br1_tr1", "het_ut_br1_tr2"),
  wt_7d = c("wt_7d_abeta_br2_tr1", "wt_7d_abeta_br2_tr2"),
  wt_ut = c("wt_ut_br2_tr1", "wt_ut_br2_tr2")
)

X_tmp <- lapply(1:(length(group_list)-1), function(k){
  tmp <- rownames(Y)
  vec <- rep(FALSE, length(tmp))
  vec[tmp %in% group_list[[k]]] <- TRUE
  
  as.numeric(vec)
})
X <- do.call(cbind, X_tmp)
rownames(X) <- rownames(Y)
colnames(X) <- names(group_list)[-length(group_list)]
X <- data.frame(X)

# https://github.com/statdivlab/radEmu/blob/main/R/emuFit.R
# I think this is not in the usual scale... I think these p-values need a very high replicates?
ch_fit <- radEmu::emuFit(Y = Y,
                         X = X,
                         compute_cis = FALSE,
                         run_score_tests = TRUE) 
ch_fit
# note: Not sure what scale the estimate is on

#################

# what about using this one?
# https://sccoda.readthedocs.io/en/latest/using_other_compositional_methods.html

