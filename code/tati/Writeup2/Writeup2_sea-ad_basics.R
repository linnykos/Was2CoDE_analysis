library(Seurat)
object <- readRDS("microglia_mtg.rds")
head(object@meta.data)
colnames(object@meta.data)

# "Cognitive status", "sex", "PMI", "donor_id", "Age at death", (no batch variable for this one)

keep_vec <- rep(TRUE, ncol(object))
keep_vec[which(object@meta.data[,"Cognitive status"] == "Reference")] <- FALSE
object$keep <- keep_vec

object <- subset(object, keep == TRUE)
table(object@meta.data[,"Cognitive status"])

table(object@meta.data[,"donor_id"], object@meta.data[,"Cognitive status"])

# when you run custom ideas, use this matrix:
mat <- SeuratObject::LayerData(object, 
                               layer = "data",
                               assay = "RNA")
# this matrix is genes-by-cells
mat[1:5,1:5]

###################

load("~/kzlinlab/projects/subject-de/out/kevin/Writeup6/Writeup6_prater_scvi-seurat.RData")

ss_data_norm
table(ss_data_norm$Pt_ID, ss_data_norm$CognitiveStatus)

# when you run custom ideas OR make the boxplots, use this matrix:
mat <- SeuratObject::LayerData(ss_data_norm, 
                               layer = "data",
                               assay = "RNA")
# this matrix is genes-by-cells
mat[1:5,1:5]

# to make a boxplot for a gene:
## make a new data frame with 3 columns (n rows, 1 row per cell)
# -- one column for: gene expression
# -- another: donor
# -- another: CognitiveStatus
# pass this into ggplot's boxplot <- color the AD by red, non-AD by blue
# all the red boxplots are on one side, blues on the other side

