rm(list=ls())
library(openxlsx)

metadata <- openxlsx::read.xlsx("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_cohort_donor_metadata_082222.xlsx")
rownames(metadata) <- metadata$Donor.ID
harmonized_scores <- openxlsx::read.xlsx("~/Dropbox/Collaboration-and-People/alzheimer/data/joey-harmonized-cognitive/sea_ad_cognitive_slopes_79.xlsx")
rownames(harmonized_scores) <- harmonized_scores$donor_name
neuropath <- read.csv("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv")
rownames(neuropath) <- neuropath$Donor.ID

path_idx <- grep("^percent.*area_Grey", colnames(neuropath))
colnames(neuropath)[path_idx]
neuropath <- neuropath[,path_idx]
summary(harmonized_scores)

included_idx <- which(metadata$Donor.ID %in% harmonized_scores$donor_name)
metadata <- metadata[included_idx,]

ad_idx <- intersect(
  which(metadata$Braak %in% c("Braak V", "Braak VI")),
  which(metadata$CERAD.score %in% c("Moderate", "Frequent"))
)
metadata <- metadata[ad_idx,]
harmonized_scores <- harmonized_scores[rownames(metadata),]
neuropath <- neuropath[rownames(metadata),]

#################

neuropath2 <- scale(neuropath)
# flip specifically flip NeuN
idx <- grep("NeuN", colnames(neuropath2))
neuropath2[,idx] <- -neuropath2[,idx]

metadata2 <- metadata[,c("Age.at.Death", "Sex", "Years.of.education", "APOE4.Status")]
metadata2[which(metadata2[,"Age.at.Death"] == "90+"),"Age.at.Death"] <- "90"
metadata2[,"Age.at.Death"] <- as.numeric(metadata2[,"Age.at.Death"])
metadata2[,"Sex"] <- as.numeric(as.factor(metadata2[,"Sex"]))
metadata2[,"APOE4.Status"] <- as.numeric(as.factor(metadata2[,"APOE4.Status"]))
df <- as.data.frame(cbind(-harmonized_scores$slope_zmem0, neuropath2, 
                          metadata2))
colnames(df)[1] <- "y"
lm_res <- stats::lm(y ~ ., data = df)
stats::coef(lm_res)

plot(x = lm_res$fitted.values,
     y = -harmonized_scores$slope_zmem0, asp = T)
lines(x = c(-10,10), y = c(-10,10), col = "red", lwd = 2, lty = 2)
