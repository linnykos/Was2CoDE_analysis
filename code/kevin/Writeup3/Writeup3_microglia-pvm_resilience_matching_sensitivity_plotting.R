rm(list=ls())
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching_sensitivity.RData")

sa_mat <- sa_mat_safe
for(i in 1:nrow(sa_mat)){
  sa_mat[i,] <- sort(sa_mat[i,],decreasing = F)
}
sa_mat2 <- sa_mat[c("XRCC2", "EFCAB6", "TNS1", "BAP1"),]

ylim <- range(as.numeric(abs(sa_mat2)))

png("~/kzlinlab/projects/subject-de/git/subject-de_kevin/figures/kevin/Writeup3/Writeup3_sensitivity-analysis.png",
    height = 3000, width = 3000, units = "px", res = 300)
plot(NA, xlim = c(1, 6), ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n")
