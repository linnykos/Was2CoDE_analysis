rm(list=ls())
n <- 20
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(n)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(n)
transparent_gray <- rgb(0.5,0.5,0.5,0.5)
two_letters <- substr(transparent_gray, start = 8, stop = 9)
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

set.seed(10)
case_mean_vec <- 0.3 + runif(n, min = -0.2, 0.2)
case_sd_vec <- rep(0.3, n)
case_size_vec <- rep(1, n)

control_mean_vec <- 0.3 + runif(n, min = -0.2, 0.2)
control_mean_vec[1] <- 1
control_sd_vec <- rep(0.3, n); # control_sd_vec[1] <- 0.1
control_size_vec <- rep(1, n); control_size_vec[1] <- 4

xlim <- c(-1,2)
ylim <- c(0, 20)

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

# png(paste0("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/git/subject-de/figures/kevin/Writeup5/de-illustration_bad-singlecell.png"),
#     height = 1200, width = 2000,
#     units = "px", res = 500)
par(mar = c(3,0.25,0,0.25), bg = NA)
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
# graphics.off()

#################

# actually test out this seting
set.seed(10)
case_list <- lapply(1:n, function(i){
  stats::rnorm(n = case_size_vec[i]*100,
               mean = case_mean_vec[i],
               sd = case_sd_vec[i])
})
control_list <- lapply(1:n, function(i){
  stats::rnorm(n = control_size_vec[i]*100,
               mean = control_mean_vec[i],
               sd = control_sd_vec[i])
})

# single-cell testing
case_cells <- unlist(case_list)
control_cells <- unlist(control_list)
stats::wilcox.test(x = case_cells,
                   y = control_cells)

# donor testing
avg_posterior_case_mean_mat <- sapply(case_list, mean)
avg_posterior_case_var_mat <- sapply(case_list, stats::var)
avg_posterior_control_mean_mat <- sapply(control_list, mean)
avg_posterior_control_var_mat <- sapply(control_list, stats::var)

case_gaussian_mean <- mean(avg_posterior_case_mean_mat)
control_gaussian_mean <- mean(avg_posterior_control_mean_mat)
case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = matrix(avg_posterior_case_mean_mat, ncol = 1),
  avg_posterior_var_mat = matrix(avg_posterior_case_var_mat, ncol = 1)
)
control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = matrix(avg_posterior_control_mean_mat, ncol = 1),
  avg_posterior_var_mat = matrix(avg_posterior_control_var_mat, ncol = 1)
)

n1 <- n
n2 <- n
teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
  (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))

df_num <- (case_gaussian_var/n1 + control_gaussian_var/n2)^2
df_denom <- (case_gaussian_var/n1)^2/(n1-1) + (control_gaussian_var/n2)^2/(n2-1)
df_vec <- df_num/df_denom

1-stats::pt(abs(teststat_vec), df = df_vec)





