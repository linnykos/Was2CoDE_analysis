load("~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideas_subset.RData")
library(dplyr)

AD_patients <- subset(meta_ind, CognitiveStatus == "Dementia")
no_AD_patients <- subset(meta_ind, CognitiveStatus == "No_dementia")
set.seed(10)
selected_AD <- AD_patients[sample(nrow(AD_patients), 2), ]
selected_no_AD <- no_AD_patients[sample(nrow(no_AD_patients), 2), ]
selected_patients <- rbind(selected_AD, selected_no_AD)
meta_ind <- selected_patients
meta_cell <- meta_cell[meta_cell$Pt_ID %in% selected_patients$individual, ]
count_matrix <- count_matrix_subset[, colnames(count_matrix_subset) %in% rownames(meta_cell)]
save(count_matrix, meta_cell, meta_ind, var_per_cell, var2test, var2test_type,
     file = "~/kzlinlab/projects/subject-de/out/tati/Writeup1/Writeup1_microglia_ideas_subset_testing.RData")
