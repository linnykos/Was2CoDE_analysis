rm(list=ls())
n <- 10
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(n)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(n)
transparent_gray <- rgb(0.5,0.5,0.5,0.3)
two_letters <- substr(transparent_gray, start = 8, stop = 9)
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

set.seed(10)
case_mean_vec <- 0.6 + runif(n, min = -0.3, 0.3)
case_sd_vec <- rep(0.3, n)
case_size_vec <- rep(1, n)

control_mean_vec <- 0.3 + runif(n, min = -0.3, 0.3)
control_sd_vec <- rep(0.3, n)
control_size_vec <- rep(1, n)

xlim <- c(-1,2)
ylim <- c(0, 10)

case_gaussian_list <- lapply(1:n, function(i){
  mean_val <- case_mean_vec[i]
  sd_val <- case_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4*case_size_vec[i]
  cbind(xseq, yseq)
})
control_gaussian_list <- lapply(1:n, function(i){
  mean_val <- control_mean_vec[i]
  sd_val <- control_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4*control_size_vec[i]
  cbind(xseq, yseq)
})

par(mar = c(3,0.25,0,0.25))
plot(NA,
     xlim = xlim,
     ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
for(i in 1:n){
  graphics::polygon(x = c(control_gaussian_list[[i]][,1], rev(control_gaussian_list[[i]][,1])),
                    y = 2+c(control_gaussian_list[[i]][,2], rep(0, nrow(control_gaussian_list[[i]]))),
                    col = control_color_trans_palette[i])
}
for(i in 1:n){
  graphics::polygon(x = c(case_gaussian_list[[i]][,1], rev(case_gaussian_list[[i]][,1])),
                    y = c(case_gaussian_list[[i]][,2], rep(0, nrow(case_gaussian_list[[i]]))),
                    col = case_color_trans_palette[i])
}

for(i in 1:n){
  graphics::rug(control_mean_vec[i],
                col = control_color_palette[i],
                lwd = 2)
  graphics::rug(case_mean_vec[i],
                col = case_color_palette[i],
                lwd = 2)
}

axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
