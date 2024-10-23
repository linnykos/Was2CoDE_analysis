rm(list=ls())
library(Seurat)

out_folder <- "~/Dropbox/Collaboration-and-People/tati/out/kevin/Writeup13/"
plot_folder <- "~/Dropbox/Collaboration-and-People/tati/git/subject-de/figures/kevin/Writeup13/"
load(paste0(out_folder, "Writeup13_esvd.RData"))

true_vec <- (true_fdr_vec <= 0.5)

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

set.seed(10)
prob_vec <- rank(eSVD_obj$pvalue_list$log10pvalue)
prob_vec <- (length(prob_vec) - prob_vec)/length(prob_vec)

shuff_percentage <- 0.1
prob_vec_lowShuffle <- prob_vec
target_idx <- sample(1:length(prob_vec), size = round(shuff_percentage * length(prob_vec)))
prob_vec_lowShuffle[target_idx] <- prob_vec_lowShuffle[sample(target_idx)]
tol <- 0.1
prob_vec_lowShuffle <- prob_vec_lowShuffle + runif(length(prob_vec), min = -tol, max = tol)
prob_vec_lowShuffle <- (prob_vec_lowShuffle - min(prob_vec_lowShuffle))/diff(range(prob_vec_lowShuffle))

shuff_percentage <- 0.5
prob_vec_highShuffle <- prob_vec
target_idx <- sample(1:length(prob_vec), size = round(shuff_percentage * length(prob_vec)))
prob_vec_highShuffle[target_idx] <- prob_vec_highShuffle[sample(target_idx)]

roc_res <- roc_func(true_vec = true_vec,
                    prob_vec = prob_vec)
roc_res_low <- roc_func(true_vec = true_vec,
                    prob_vec = prob_vec_lowShuffle)
roc_res_high <- roc_func(true_vec = true_vec,
                    prob_vec = prob_vec_highShuffle)

png(paste0(plot_folder, "Writeup3_simulation-power.png"),
    height = 1000, width = 800, res = 300, units = "px")
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
lines(y = roc_res_high[,2],
      x = roc_res_high[,3],
      lwd = 4,
      col = "black")
lines(y = roc_res[,2],
      x = roc_res[,3],
      lwd = 2.5,
      lty = 2,
      col = "black")
lines(y = roc_res_low[,2],
      x = roc_res_low[,3],
      lwd = 4,
      col = rgb(181, 198, 5, maxColorValue = 255))
dev.off()
