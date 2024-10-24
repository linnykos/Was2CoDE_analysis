rm(list=ls())
library(foreach)
library(future)
library(rngtools)
library(Seurat)
set.seed(10)

################################################
############Stratifying data for plotting ###############
################################################

############ Katie's data ###############
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_prater_scVI-postprocessed.RData")
summary(ss_data_norm@meta.data)

# average read depth
median(ss_data_norm$nCount_RNA)

# average number of cells per donor
median(table(ss_data_norm$Pt_ID))

metadata <- ss_data_norm@meta.data
metadata_subset <- metadata[,c("Pt_ID", "Study_Designation", "Sex", "APOEe4_status","Race","PMI","coded_Age", "SeqBatch")]
metadata_subset <- unique(metadata_subset)
dim(metadata_subset) # one row per unique donor
table(metadata_subset$Study_Designation, metadata_subset$Sex) 
table(metadata_subset$Study_Designation, metadata_subset$APOEe4_status)
table(metadata_subset$Study_Designation, metadata_subset$coded_Age) 

############ SEA-AD###############
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup10/Writeup10_sea-ad_microglia_scVI-postprocessed.RData") 
summary(seurat_obj@meta.data)

# average read depth
median(seurat_obj$nCount_RNA)

# average number of cells per donor
median(table(seurat_obj$donor_id))

metadata <- seurat_obj@meta.data
metadata_subset <- metadata[,c("donor_id", "ADNC", "sex", "APOE4status","self_reported_ethnicity","PMI","Ageatdeath")]
metadata_subset <- unique(metadata_subset)
dim(metadata_subset) # one row per unique donor
table(metadata_subset$ADNC, metadata_subset$sex) 
table(metadata_subset$ADNC, metadata_subset$APOE4status)
table(metadata_subset$ADNC, metadata_subset$Ageatdeath) 

############ ROSMAP ###############
load("~/kzlinlab/projects/subject-de/out/kevin/Writeup11/Writeup11_rosmap_scVI-postprocessed.RData") 
summary(seurat_obj@meta.data)

# average read depth
median(seurat_obj$nCount_RNA)

# average number of cells per donor
median(table(seurat_obj$Pt_ID))

metadata <- seurat_obj@meta.data
metadata_subset <- metadata[,c("Pt_ID", "ADpath", "sex", "APOEe4_status","race","age_death")]
metadata_subset <- unique(metadata_subset)
dim(metadata_subset) # one row per unique donor
table(metadata_subset$ADpath, metadata_subset$sex) 
table(metadata_subset$ADpath, metadata_subset$APOEe4_status)
table(metadata_subset$ADpath, metadata_subset$age_death)

metadata_subset$age_group <- cut(
  metadata_subset$age_death,
  breaks = c(-Inf, 70, 85, Inf),  # Define the age brackets: <=70, 71-85, >=86
  labels = c("<=70", "71-85", ">=86"),  
  right = TRUE  #  intervals includes upper bound
)

age_adpath_counts <- table(metadata_subset$ADpath, metadata_subset$age_group)
age_adpath_counts

####################################
############ Plotting###############
####################################
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(reshape2)
library(scales) 

############ Prater's ###############

# Data for gender distributio
gender_data <- data.frame(
  Group = c('Case', 'Control'),
  Female = c(9, 6),
  Male = c(3, 4)
)

# Data for APOE4 distribution
apoe_data <- data.frame(
  Group = c('Case', 'Control'),
  No = c(6, 9),
  Yes = c(6, 1)
)

age_data <- data.frame(
  Group = c('Case', 'Control'),
  `>=86` = c(8, 5),      
  `71-85` = c(3, 4),
  `<=70` = c(1, 1),
  check.names = FALSE  # Disable automatic conversion of column names
)

create_stacked_bar <- function(data, title, fill_label) {
  data_long <- reshape2::melt(data, id.vars = 'Group')
  data_long$variable <- factor(data_long$variable)
  ggplot(data_long, aes(y = Group, x = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
    labs(x = "Number of Donors", 
         y = NULL, 
         fill = fill_label,
         title = title) +
    scale_x_continuous(
      limits = c(0, 12),
      breaks = seq(0, 12, by = 2),
      expand = c(0, 0) 
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
    )
}

gender_plot <- create_stacked_bar(gender_data, "Gender Distribution", "Gender") +
  scale_fill_manual(values = c("Female" = "#FF9999", "Male" = "#66B2FF"))

apoe_plot <- create_stacked_bar(apoe_data, "APOE4 Status Distribution", "APOE4") +
  scale_fill_manual(values = c("No" = "#99FF99", "Yes" = "#FF99FF"))

age_plot <- create_stacked_bar(age_data, "Age Distribution", "Age Group") +
  scale_fill_manual(values = c(
    '>=86' = "#FCC351",  
    '71-85' = "#A1DEE0", 
    '<=70' = "#FD8D6E" 
  ))
####################################

############ SEA-AD ###############
# Data for gender distribution
gender_data <- data.frame(
  Group = c('Case', 'Control'),
  Female = c(36, 12),
  Male = c(23, 9) 
)

# Data for APOE4 distribution
apoe_data <- data.frame(
  Group = c('Case', 'Control'),
  No = c(36, 21),
  Yes = c(21, 0)
)

age_data <- data.frame(
  Group = c('Case', 'Control'),
  `>=86` = c(45, 12),      
  `71-85` = c(11, 9),
  `<=70` = c(3, 0),
  check.names = FALSE  # Disable automatic conversion of column names
)

create_stacked_bar <- function(data, title, fill_label) {
  data_long <- reshape2::melt(data, id.vars = 'Group')
  data_long$variable <- factor(data_long$variable)
  ggplot(data_long, aes(y = Group, x = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
    labs(x = "Number of Donors", 
         y = NULL, 
         fill = fill_label,
         title = title) +
    scale_x_continuous(
      limits = c(0, 59),
      breaks = seq(0, 59, by = 10),
      expand = c(0, 0) 
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      )
}

gender_plot <- create_stacked_bar(gender_data, "Gender Distribution", "Gender") +
  scale_fill_manual(values = c("Female" = "#FF9999", "Male" = "#66B2FF"))

apoe_plot <- create_stacked_bar(apoe_data, "APOE4 Status Distribution", "APOE4") +
  scale_fill_manual(values = c("No" = "#99FF99", "Yes" = "#FF99FF"))

age_plot <- create_stacked_bar(age_data, "Age Distribution", "Age Group") +
  scale_fill_manual(values = c(
    '>=86' = "#FCC351",  
    '71-85' = "#A1DEE0", 
    '<=70' = "#FD8D6E" 
  ))
################################

############ ROSMAP ###############
# Data for gender distribution 
gender_data <- data.frame(
  Group = c('Case', 'Control'),
  Female = c(101, 74),
  Male = c(91, 79)  
)

# Data for APOE4 distribution
apoe_data <- data.frame(
  Group = c('Case', 'Control'),
  No = c(139, 116),
  Yes = c(53, 37)
)

age_data <- data.frame(
  Group = c('Case', 'Control'),
  `>=86` = c(136, 104),      
  `71-85` = c(56, 49),
  `<=70` = c(0, 0),
  check.names = FALSE  # Disable automatic conversion of column names
)

age_data2 <- data.frame(
  Group = c('Case', 'Control'),
  `>=86` = c(136, 104),      
  `71-85` = c(56, 49),
  `<=70` = c(0, 0),
  check.names = FALSE  # Disable automatic conversion of column names
)

create_stacked_bar <- function(data, title, fill_label) {
  data_long <- reshape2::melt(data, id.vars = 'Group')
  data_long$variable <- factor(data_long$variable)
  ggplot(data_long, aes(y = Group, x = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
    labs(x = "Number of Donors", 
         y = NULL, 
         fill = fill_label,
         title = title) +
    scale_x_continuous(
      limits = c(0, 192),
      breaks = seq(0, 192, by = 25),
      expand = c(0, 0) 
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
    )
}

gender_plot <- create_stacked_bar(gender_data, "Gender Distribution", "Gender") +
  scale_fill_manual(values = c("Female" = "#FF9999", "Male" = "#66B2FF"))

apoe_plot <- create_stacked_bar(apoe_data, "APOE4 Status Distribution", "APOE4") +
  scale_fill_manual(values = c("No" = "#99FF99", "Yes" = "#FF99FF"))

age_plot <- create_stacked_bar(age_data, "Age Distribution", "Age Group") +
  scale_fill_manual(values = c(
    '>=86' = "#FCC351",  
    '71-85' = "#A1DEE0", 
    '<=70' = "#FD8D6E" 
  ))
#############################################
#############################################
#############################################
all_datasets <- data.frame(
  Dataset = rep(c("Prater", "SEA-AD", "ROSMAP"), each = 4),
  Metric = rep(c("Total Cells", "Median Cells per Donor", "Average Read Depth", "Number of Donors"), 3),
  Value = c(
    # Prater's 
    127367, 5764.5, 584, 22,
    # SEA-AD 
    41187, 496, 1476, 80,
    # ROSMAP 
    71015, 186, 2463, 345
  )
)

all_datasets$Metric <- factor(
  all_datasets$Metric, 
  levels = c("Number of Donors", "Average Read Depth", "Median Cells per Donor", "Total Cells")  # Set order
)

all_datasets <- all_datasets %>%
  mutate(
    Label = ifelse(Metric %in% c("Total Cells", "Number of Donors"), 
                   format(as.integer(Value), big.mark = ",", scientific = FALSE), 
                   format(Value, big.mark = ",", scientific = FALSE))
  )

# Calculate normalized values within each dataset relative to the maximum value in that metric
all_datasets <- all_datasets %>%
  group_by(Metric) %>%
  mutate(NormalizedValue = Value / max(Value)) %>%
  ungroup()

create_dataset_plot <- function(data, dataset_name) {
  dataset_data <- data %>% filter(Dataset == dataset_name)
  
  ggplot(dataset_data, aes(y = Metric, x = NormalizedValue, fill = Metric)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = Label), 
              hjust = -0.05,
              size = 2) +
    labs(title = paste(dataset_name, "Dataset"),
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 70, b = 20, l = 20)
    ) +
    scale_x_continuous(
      limits = c(0, 1.2),  # Reduced upper limit to make bars longer
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = c(
      "Total Cells" = "#fbf1d1",          
      "Median Cells per Donor" = "#fcdfe5", 
      "Average Read Depth" = "#daf1ee",   
      "Number of Donors" = "#b6e3e7"
    ))
}

prater_plot <- create_dataset_plot(all_datasets, "Prater")
seaad_plot <- create_dataset_plot(all_datasets, "SEA-AD")
rosmap_plot <- create_dataset_plot(all_datasets, "ROSMAP")