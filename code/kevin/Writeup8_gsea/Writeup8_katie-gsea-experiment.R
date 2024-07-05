# from /Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/git/fxomics/13_snRNAseq_diff_gene_expression_v5_clust.Rmd

KEGG_pathways <- read.gmt("../data/c2.cp.kegg.v7.2.symbols.gmt")
React_pathways <- read.gmt("../data/c2.cp.reactome.v7.2.symbols.gmt")
WP_pathways <- read.gmt("../data/c2.cp.wikipathways.v7.2.symbols.gmt")
GO_pathways <- read.gmt("../data/c5.go.v7.2.symbols.gmt")

pathway_list <- c("KEGG_pathways",
                  "React_pathways",
                  "WP_pathways",
                  "GO_pathways")

group_markers <- as.data.table(group_markers, keep.rownames = TRUE) %>%
  dplyr::select(gene, avg_logFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene) %>%
  deframe() %>%
  sort(decreasing = T)

# Run the gsea algorithm for each pathway
for (pathway in pathway_list) {
  gsea_results <- run_gsea(get(pathway), group_markers)
  if (dim(gsea_results@result)[1] > 0) {
    save_gsea(
      paste0(cluster_idents[comp_clust[cluster] +
                              1 - as.integer(leiden_alg)], "_vs_", cluster_idents[base_clust[index] + 1  -
                                                                                    as.integer(leiden_alg)]),
      gsea_results,
      pathway,
      chose_resolution
    )
  }
}