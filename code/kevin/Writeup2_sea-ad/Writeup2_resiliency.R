rm(list=ls())
library(openxlsx)

metadata <- openxlsx::read.xlsx("~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/sea-ad_cohort_donor_metadata_082222.xlsx")

donor_id <- metadata$Donor.ID
write.csv(donor_id, 
          file = "~/Dropbox/Collaboration-and-People/alzheimer/data/sea-ad/donor_id.csv")

# extract the specific columns I need
metadata_ss <- metadata[,c("Last.CASI.Score",
                           "Interval.from.last.CASI.in.months",
                           "Thal",
                           "Braak",
                           "CERAD.score",
                           "Overall.CAA.Score",
                           "Total.microinfarcts.in.screening.sections",
                           "LATE")]

# preprocess certain columns
## Braak
tmp <- factor(metadata_ss$Braak, levels = c("Braak 0", "Braak II", "Braak III", "Braak IV", "Braak V", "Braak VI"))
tmp2 <- as.numeric(tmp)
metadata_ss$Braak <- tmp2

table(metadata$Braak, metadata$CERAD.score)
table(metadata$Braak, metadata$Cognitive.Status)
table(metadata$Braak, metadata$Last.CASI.Score)

idx <- intersect(
  intersect(which(metadata$Braak %in% c("Braak V", "Braak VI")),
                 which(metadata$CERAD.score %in% c("Moderate", "Frequent"))),
  which(metadata$Interval.from.last.CASI.in.months <= 36)
)
length(idx)
hist(metadata$Last.CASI.Score[idx])

## Thal
tmp <- factor(metadata_ss$Thal, levels = c("Thal 0", "Thal 1", "Thal 2", "Thal 3", "Thal 4", "Thal 5"))
tmp2 <- as.numeric(tmp)
metadata_ss$Thal <- tmp2

## CERAD.score
tmp <- factor(metadata_ss$CERAD.score, levels = c("Absent", "Sparse", "Moderate", "Frequent"))
tmp2 <- as.numeric(tmp)
metadata_ss$CERAD.score <- tmp2

## Overall.CAA.Score
tmp <- factor(metadata_ss$Overall.CAA.Score, levels = c("Not identified", "Mild", "Moderate", "Severe"))
tmp2 <- as.numeric(tmp)
metadata_ss$Overall.CAA.Score <- tmp2

## CERAD.score
tmp <- factor(metadata_ss$CERAD.score, levels = c("1", "2", "3", "4"))
tmp2 <- as.numeric(tmp)
metadata_ss$CERAD.score <- tmp2

## LATE
tmp <- factor(metadata_ss$LATE, levels = c("Not identified", "LATE Stage 1", "LATE Stage 2", "LATE Stage 3", "Unclassifiable"))
tmp2 <- as.numeric(tmp)
metadata_ss$LATE <- tmp2

# remove all patients that had their last score more than 24 months before death
idx <- which(metadata_ss[,"Interval.from.last.CASI.in.months"] <= 24)
metadata_ss <- metadata_ss[idx,]

# replace all the NA's with the mean
for(j in 1:ncol(metadata_ss)){
  if(any(is.na(metadata_ss[,j]))){
    idx <- which(is.na(metadata_ss[,j]))
    metadata_ss[idx,j] <- mean(metadata_ss[,j], na.rm = T)
  }
}

# now do a nonnegative least squres
lm_res <- stats::lm(
  Last.CASI.Score ~ Thal + Braak + CERAD.score + Overall.CAA.Score + Total.microinfarcts.in.screening.sections + LATE,
  data = metadata_ss
)

plot(x = lm_res$fitted.values,
     y = metadata_ss[,"Last.CASI.Score"],
     asp = T, pch = 16)
lines(x = c(-1e4,1e4), y = c(-1e4,1e4), col = 2, lwd = 2, lty = 2)

summary(lm_res)
