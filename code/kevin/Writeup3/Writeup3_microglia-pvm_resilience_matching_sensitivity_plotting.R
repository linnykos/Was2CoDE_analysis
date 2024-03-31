rm(list=ls())
load("../../../out/kevin/Writeup3/Writeup3_sea-ad_microglia_matching_sensitivity.RData")

sa_mat <- sa_mat_safe
for(i in 1:nrow(sa_mat)){
  sa_mat[i,] <- sort(sa_mat[i,],decreasing = F)
}
sa_mat2 <- abs(sa_mat[c("XRCC2", "EFCAB6", "TNS1", "BAP1"),])

ylim <- range(as.numeric(abs(sa_mat2)))

orange_col <- rgb(235, 134, 47, maxColorValue = 255) 
blue_col <- rgb(50, 174, 255, maxColorValue = 255) 
green_col <- rgb(70, 177, 70, maxColorValue = 255) 

col_vec <- c("black", orange_col, blue_col, green_col)

png("../figures/kevin/Writeup3/Writeup3_sensitivity-analysis.png",
   height = 1200, width = 1200, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot(NA, 
     bty = "n",
     xaxt = "n", 
     xlab = "",
     xlim = c(1, 6), 
     yaxt = "n", 
     ylab = "",
     ylim = ylim)
axis(1); axis(2)
y_breaks <- seq(0.1, 0.5, by = 0.05)
for(y in y_breaks){
  lines(x = c(-10,10), y = rep(y,2), col = "gray", lty = 2)
}
for(i in 1:4){
  points(x = 1:6, y = sa_mat2[i,], pch = 16, col = col_vec[i], cex = 2.5)
  lines(x = 1:6, y = sa_mat2[i,], lwd = 2, col = col_vec[i])
}
graphics.off()


