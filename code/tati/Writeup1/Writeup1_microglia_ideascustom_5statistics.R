# Define the compute_statistics function
compute_statistics <- function(dist_list, meta_ind, pval_list) {
  n_gene <- dim(dist_list$distance)[1]
  n_donors <- dim(dist_list$distance)[2]
  
  statistics <- matrix(NA, nrow = n_gene, ncol = 5)
  colnames(statistics) <- c("median_loc_dist2_case_control", 
                            "median_loc_dist2_case_case_control_control", 
                            "median_size_dist2_case_control", 
                            "median_size_dist2_case_case_control_control", 
                            "-log10_p_value")
  
  for (i in 1:n_gene) {
    location <- dist_list$location[i, , ]
    size <- dist_list$size[i, , ]
    distance <- dist_list$distance[i, , ]
    
    loc_dist2 <- location / (distance^2)
    size_dist2 <- size / (distance^2)
    
    case_indices <- which(meta_ind$CognitiveStatus == "Dementia")
    control_indices <- which(meta_ind$CognitiveStatus == "No_dementia")
    
    loc_dist2_case_control <- loc_dist2[case_indices, control_indices]
    size_dist2_case_control <- size_dist2[case_indices, control_indices]
    
    loc_dist2_case_case_control_control <- loc_dist2[c(case_indices, control_indices), c(case_indices, control_indices)]
    size_dist2_case_case_control_control <- size_dist2[c(case_indices, control_indices), c(case_indices, control_indices)]
    
    statistics[i, 1] <- median(loc_dist2_case_control, na.rm = TRUE)
    statistics[i, 2] <- median(loc_dist2_case_case_control_control, na.rm = TRUE)
    statistics[i, 3] <- median(size_dist2_case_control, na.rm = TRUE)
    statistics[i, 4] <- median(size_dist2_case_case_control_control, na.rm = TRUE)
    
    # Assuming we want to use the p-value for the 'distance' metric for the log transformation
    p_value <- pval_list[[1]][i]
    statistics[i, 5] <- -log10(p_value)
  }
  
  return(statistics)
}

# Unit test
library(testthat)
# Mock data to simulate the dist_list, meta_ind, and pval_list
mock_dist_list <- list(
  distance = array(runif(2 * 4 * 4), dim = c(2, 4, 4)),
  location = array(runif(2 * 4 * 4), dim = c(2, 4, 4)),
  location_sign = array(runif(2 * 4 * 4), dim = c(2, 4, 4)),
  size = array(runif(2 * 4 * 4), dim = c(2, 4, 4)),
  size_sign = array(runif(2 * 4 * 4), dim = c(2, 4, 4)),
  shape = array(runif(2 * 4 * 4), dim = c(2, 4, 4))
)

mock_meta_ind <- data.frame(
  individual = c("ind1", "ind2", "ind3", "ind4"),
  CognitiveStatus = c("Dementia", "No_dementia", "Dementia", "No_dementia")
)

mock_pval_list <- list(
  runif(2),
  runif(2),
  runif(2),
  runif(2),
  runif(2),
  runif(2)
)

test_that("compute_statistics works correctly", {
  result <- compute_statistics(mock_dist_list, mock_meta_ind, mock_pval_list)
  
  expected_result <- matrix(NA, nrow = 2, ncol = 5)
  expected_result[1, 1] <- median(mock_dist_list$location[1, c(1, 3), c(2, 4)] / mock_dist_list$distance[1, c(1, 3), c(2, 4)]^2, na.rm = TRUE)
  expected_result[1, 2] <- median(mock_dist_list$location[1, , ] / mock_dist_list$distance[1, , ]^2, na.rm = TRUE)
  expected_result[1, 3] <- median(mock_dist_list$size[1, c(1, 3), c(2, 4)] / mock_dist_list$distance[1, c(1, 3), c(2, 4)]^2, na.rm = TRUE)
  expected_result[1, 4] <- median(mock_dist_list$size[1, , ] / mock_dist_list$distance[1, , ]^2, na.rm = TRUE)
  expected_result[1, 5] <- -log10(mock_pval_list[[1]][1])
  
  expected_result[2, 1] <- median(mock_dist_list$location[2, c(1, 3), c(2, 4)] / mock_dist_list$distance[2, c(1, 3), c(2, 4)]^2, na.rm = TRUE)
  expected_result[2, 2] <- median(mock_dist_list$location[2, , ] / mock_dist_list$distance[2, , ]^2, na.rm = TRUE)
  expected_result[2, 3] <- median(mock_dist_list$size[2, c(1, 3), c(2, 4)] / mock_dist_list$distance[2, c(1, 3), c(2, 4)]^2, na.rm = TRUE)
  expected_result[2, 4] <- median(mock_dist_list$size[2, , ] / mock_dist_list$distance[2, , ]^2, na.rm = TRUE)
  expected_result[2, 5] <- -log10(mock_pval_list[[1]][2])
  
  colnames(expected_result) <- c("median_loc_dist2_case_control", 
                                 "median_loc_dist2_case_case_control_control", 
                                 "median_size_dist2_case_control", 
                                 "median_size_dist2_case_case_control_control", 
                                 "-log10_p_value")
  
  expect_equal(result, expected_result)
})

print("All tests passed.")

#run data
load("~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideascustom.RData")
statistics <- compute_statistics(dist_list, meta_ind, pval_ideas)
save(statistics,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideascustom_5statistics.RData")

print("Done! :)")

#Plots

library(ggplot2)
library(ggrepel)
library(dplyr)
# Convert the statistics matrix to a data frame for ggplot2
statistics_df <- as.data.frame(statistics)
# Arrange the data frame by -log10_p_value in descending order
statistics_df <- statistics_df %>%
  arrange(desc("-log10_p_value"))

colnames(statistics_df) <- c("median_loc_dist2_case_control", 
                             "median_loc_dist2_case_case_control_control", 
                             "median_size_dist2_case_control", 
                             "median_size_dist2_case_case_control_control", 
                             "-log10_p_value")

# Labels
statistics_df$label <- ifelse(rank(-statistics_df$"-log10_p_value") <= 10, rownames(statistics_df), "")

# Plot 1: Median of location/distance^2 for pairs of donors case-control vs
# Median of location/distance^2 for pairs of donors case-case or control-control
plot1 <- ggplot(statistics_df, aes(x = median_loc_dist2_case_case_control_control, 
                                   y = median_loc_dist2_case_control, 
                                   label = label, 
                                   color = `-log10_p_value`)) +
  geom_point() +
  geom_text_repel(max.overlaps = 200) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 1: Med location/distance^2 (Case-Control vs Case-Case or Control-Control)",
       x = "Med location/dist^2 (Case-Case or Control-Control)",
       y = "Med location/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  scale_color_gradient(low = "#A9C8F3", high = "#0C2389")

#print(plot1)
ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de_location_distance.png"),
                plot1, device = "png", width = 5, height = 7, units = "in")

# Plot 2: Median of size/distance^2 for pairs of donors case-control vs
# Median of size/distance^2 for pairs of donors case-case or control-control
plot2 <- ggplot(statistics_df, aes(x = median_size_dist2_case_case_control_control, 
                                   y = median_size_dist2_case_control, 
                                   label = label, 
                                   color = `-log10_p_value`)) +
  geom_point() +
  geom_text_repel(max.overlaps = 200) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 2: Median of size/distance^2 (Case-Control vs Case-Case or Control-Control)",
       x = "Median size/dist^2 (Case-Case or Control-Control)",
       y = "Median size/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_color_gradient(low = "#A9C8F3", high = "#0C2389")

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de_size_distance.png"),
                plot2, device = "png", width = 5, height = 7, units = "in")


# Plot 3: Median of size/distance^2 for pairs of donors case-control vs
# Median of location/distance^2 for pairs of donors case-control
plot3 <- ggplot(statistics_df, aes(x = median_loc_dist2_case_control, 
                                   y = median_size_dist2_case_control, 
                                   label = label,
                                   color = `-log10_p_value`)) +
  geom_point() +
  geom_text_repel(max.overlaps = 200) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Plot 3: size/distance^2 vs location/dist^2 Case-Control",
       x = "Median size/dist^2 (Case-Control)",
       y = "Median location/dist^2 (Case-Control)",
       color = "-log10(p-value)") +
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_color_gradient(low = "#A9C8F3", high = "#0C2389")

ggplot2::ggsave(filename = paste0("~/kzlinlab/projects/subject-de/git/subject-de_tati/figures/tati/Writeup1/Writeup1_de_size_location.png"),
                plot3, device = "png", width = 5, height = 7, units = "in")
