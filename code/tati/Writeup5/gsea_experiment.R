load("~/kzlinlab/projects/subject-de/out/tati/Writeup5/Writeup5_prater_NEBULA.RData")
head(nebula_res$summary)

teststat_vec <- -log10(nebula_res$summary[,"p_Study_DesignationAD"]) # our "proxy" for a test-statistics
names(teststat_vec) <- nebula_res$summary[,"gene"]

# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("AnnotationDbi")

library(org.Hs.eg.db)
library(clusterProfiler)

gse <- clusterProfiler::gseGO(
  teststat_vec,
  ont = "BP", # what kind of pathways are you interested in
  keytype = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,       # p-value threshold for pathways
  minGSSize = 10,            # minimum gene set size
  maxGSSize = 500            # maximum gene set size
)