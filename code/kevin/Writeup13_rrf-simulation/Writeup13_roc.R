rm(list=ls())
library(Seurat)

out_folder <- "~/Dropbox/Collaboration-and-People/tati/out/kevin/Writeup13/"
load(paste0(out_folder, "Writeup13_esvd.RData"))

true_vec <- (true_fdr_vec <= 0.05)

# compute the ROC curves
roc_func <- function(true_vec, prob_vec){
  stopifnot(all(prob_vec >= 0), all(prob_vec <= 1),
            all(true_vec %in% c(0,1)))
  seq_val <- seq(0, 1, length.out = 100)
  true_1_idx <- which(true_vec == 1)
  true_0_idx <- which(true_vec == 0)
  
  tpr_vec <- sapply(seq_val, function(x){
    bool_vec <- (prob_vec <= x)
    length(intersect(true_1_idx, which(bool_vec == 1)))/length(true_1_idx)
  })
  fpr_vec <- sapply(seq_val, function(x){
    bool_vec <- (prob_vec <= x)
    length(intersect(true_0_idx, which(bool_vec == 1)))/length(true_0_idx)
  })
  
  cbind(seq_val, tpr_vec, fpr_vec)
}

prob_vec <- rank(eSVD_obj$pvalue_list$log10pvalue)
prob_vec <- (length(prob_vec) - prob_vec)/length(prob_vec)

roc_res <- roc_func(true_vec = true_vec,
                    prob_vec = prob_vec)

plot(NA, 
     xlim = c(0,1), 
     ylim = c(0,1),
     xlab = "False positive rate",
     ylab = "True positive rate",
     xaxt = "n",
     yaxt = "n",
     bty = "n",
     asp = TRUE)
axis(1); axis(2)
lines(x = c(0,1), y = c(0,1), lty = 2)
lines(y = roc_res[,3],
      x = roc_res[,2],
      lwd = 2.5)