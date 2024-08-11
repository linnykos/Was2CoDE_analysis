rm(list=ls())

results_df <- read.table("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/git/subject-de/csv/kevin/Writeup9/Writeup9_microglia_synchrony-markers_positive_GO-analysis.txt",
                         skip = 6, 
                         sep = "\t")
colnames(results_df) <- results_df[1,]
results_df <- results_df[-1,]

# convert into numbers
colname_vec <- setdiff(colnames(results_df), c("GO biological process complete", "upload_1 (over/under)"))
for(variable in colname_vec){
  tmp <- results_df[,variable]
  tmp <- gsub(pattern = "<", 
              replacement = "", 
              x = tmp)
  results_df[,variable] <- as.numeric(tmp)
}

# remove all rows with less than 100 genes
idx <- which(results_df[,"upload_1 (951)"] >= 100)
results_df2 <- results_df[idx,]
idx <- which(results_df2[,"upload_1 (FDR)"] <= 0.05)
results_df2 <- results_df2[idx,]

# sort by pvalue
results_df3 <- results_df2[order(results_df2[,"upload_1 (FDR)"], decreasing = FALSE),]
results_df3$logpvalue <- -log10(results_df3[,"upload_1 (raw P-value)"])

results_df2 <- results_df2[order(results_df2[,"upload_1 (FDR)"], decreasing = TRUE),]

# make a plot
logpvalue <- -log10(results_df2[,"upload_1 (raw P-value)"])
num_col <- 10
break_val <- seq(min(logpvalue), max(logpvalue), length.out = num_col)
col_palette <- grDevices::colorRampPalette(c("tan", "firebrick"))(num_col)
col_vec <- sapply(logpvalue, function(x){
  col_palette[which.min(abs(x - break_val))]
})

x_vec <- results_df2[,"upload_1 (951)"]
y_vec <- results_df2[,"upload_1 (fold Enrichment)"]

plot(x = x_vec,
     y = y_vec,
     col = col_vec,
     pch = 16,
     xlab = "# genes",
     ylab = "fold enrichment")

results_df3[1:20,c("GO biological process complete", "upload_1 (951)", "upload_1 (fold Enrichment)")]

###############################
###############################
###############################

results_df <- read.table("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/tati/git/subject-de/csv/kevin/Writeup9/Writeup9_microglia_synchrony-markers_negative_GO-analysis.txt",
                         skip = 6, 
                         sep = "\t")
colnames(results_df) <- results_df[1,]
results_df <- results_df[-1,]

# convert into numbers
colname_vec <- setdiff(colnames(results_df), c("GO biological process complete", "upload_1 (over/under)"))
for(variable in colname_vec){
  tmp <- results_df[,variable]
  tmp <- gsub(pattern = "<", 
              replacement = "", 
              x = tmp)
  results_df[,variable] <- as.numeric(tmp)
}

# remove all rows with less than 100 genes
idx <- which(results_df[,"upload_1 (659)"] >= 100)
results_df2 <- results_df[idx,]
idx <- which(results_df2[,"upload_1 (FDR)"] <= 0.05)
results_df2 <- results_df2[idx,]

# sort by pvalue
results_df3 <- results_df2[order(results_df2[,"upload_1 (FDR)"], decreasing = FALSE),]
results_df3$logpvalue <- -log10(results_df3[,"upload_1 (raw P-value)"])

results_df2 <- results_df2[order(results_df2[,"upload_1 (FDR)"], decreasing = TRUE),]

# make a plot
logpvalue <- -log10(results_df2[,"upload_1 (raw P-value)"])
num_col <- 10
break_val <- seq(min(logpvalue), max(logpvalue), length.out = num_col)
col_palette <- grDevices::colorRampPalette(c("tan", "firebrick"))(num_col)
col_vec <- sapply(logpvalue, function(x){
  col_palette[which.min(abs(x - break_val))]
})

x_vec <- results_df2[,"upload_1 (659)"]
y_vec <- results_df2[,"upload_1 (fold Enrichment)"]

plot(x = x_vec,
     y = y_vec,
     col = col_vec,
     pch = 16,
     xlab = "# genes",
     ylab = "fold enrichment")

results_df3[1:20,c("GO biological process complete", 
                   "upload_1 (659)", 
                   "upload_1 (fold Enrichment)")]



