# source("Writeup5_glial_palette.R")

supertype_vec <- seurat_all$Supertype
supertype_vec[supertype_vec == "Astro_6-SEAAD"] <- "Astro_6"
supertype_vec[supertype_vec == "Micro-PVM_1_1-SEAAD"] <- "Micro-PVM_1"
supertype_vec[supertype_vec == "Micro-PVM_2_1-SEAAD"] <- "Micro-PVM_2"
supertype_vec[supertype_vec == "Micro-PVM_2_2-SEAAD"] <- "Micro-PVM_2"
supertype_vec[supertype_vec == "Micro-PVM_2_3-SEAAD"] <- "Micro-PVM_2"
supertype_vec[supertype_vec == "Micro-PVM_3-SEAAD"] <- "Micro-PVM_3"
supertype_vec[supertype_vec == "Micro-PVM_4-SEAAD"] <- "Micro-PVM_4"
supertype_vec[supertype_vec == "Oligo_2_1-SEAAD"] <- "Oligo_2"
supertype_vec[supertype_vec == "OPC_2_1-SEAAD"] <- "OPC_2"
supertype_vec[supertype_vec == "OPC_2_2-SEAAD"] <- "OPC_2"
seurat_all$Supertype2 <- supertype_vec

table(seurat_all$Supertype2)