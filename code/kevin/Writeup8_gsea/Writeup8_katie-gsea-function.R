run_gsea <- function(pathway, results) {
  user_GSEA(results, TERM2GENE = pathway, seed = T) %>%
    arrange(p.adjust) %>% filter(p.adjust < 0.05)
}

user_GSEA <- function (geneList,
                       exponent = 1,
                       minGSSize = 10,
                       maxGSSize = 500,
                       eps = 0,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       TERM2GENE,
                       TERM2NAME = NA,
                       verbose = TRUE,
                       seed = TRUE,
                       by = "fgsea",
                       ...) {
  USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
  user_GSEA_internal(
    geneList = geneList,
    exponent = exponent,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    eps = eps,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    verbose = verbose,
    USER_DATA = USER_DATA,
    seed = seed,
    by = by,
    ...
  )
}

user_GSEA_internal <- function (geneList,
                                exponent,
                                minGSSize,
                                maxGSSize,
                                eps,
                                pvalueCutoff,
                                pAdjustMethod,
                                verbose,
                                seed,
                                USER_DATA,
                                ...) {
  res <- user_fgsea(
    geneList = geneList,
    exponent = exponent,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    eps = eps,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    verbose = verbose,
    seed = TRUE,
    USER_DATA = USER_DATA,
    ...
  )
  res@organism <- "UNKNOWN"
  res@setType <- "UNKNOWN"
  res@keytype <- "UNKNOWN"
  return(res)
}

user_fgsea <- function (geneList,
                        exponent,
                        minGSSize,
                        maxGSSize,
                        eps,
                        pvalueCutoff,
                        pAdjustMethod,
                        verbose,
                        seed = TRUE,
                        USER_DATA)
{
  if (verbose)
    message("preparing geneSet collections...")
  geneSets <- get("PATHID2EXTID", envir = USER_DATA)
  if (all(!names(geneList) %in% unique(unlist(geneSets)))) {
    sg <- unlist(geneSets[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, collapse = ","))
    stop("--> No gene can be mapped....")
  }
  if (verbose)
    message("GSEA analysis...")
  tmp_res <- fgsea(
    pathways = geneSets,
    stats = geneList,
    minSize = minGSSize,
    maxSize = maxGSSize,
    eps = eps,
    gseaParam = exponent,
    nproc = 0
  )
  
  p.adj <- p.adjust(tmp_res$pval, method = pAdjustMethod)
  qvalues <- calculate_qvalue(tmp_res$pval)
  Description <- TERM2NAME(tmp_res$pathway, USER_DATA)
  params <- list(
    pvalueCutoff = pvalueCutoff,
    eps = eps,
    pAdjustMethod = pAdjustMethod,
    exponent = exponent,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  
  res <- data.frame(
    ID = as.character(tmp_res$pathway),
    Description = unname(Description),
    setSize = tmp_res$size,
    enrichmentScore = tmp_res$ES,
    NES = tmp_res$NES,
    pvalue = tmp_res$pval,
    p.adjust = p.adj,
    qvalues = qvalues,
    stringsAsFactors = FALSE
  )
  res <- res[!is.na(res$pvalue), ]
  res <- res[res$pvalue <= pvalueCutoff, ]
  res <- res[res$p.adjust <= pvalueCutoff, ]
  idx <- order(res$pvalue, decreasing = FALSE)
  res <- res[idx, ]
  if (nrow(res) == 0) {
    message("no term enriched under specific pvalueCutoff...")
    return(
      new(
        "gseaResult",
        result = res,
        geneSets = geneSets,
        geneList = geneList,
        params = params,
        readable = FALSE
      )
    )
  }
  row.names(res) <- res$ID
  observed_info <- lapply(geneSets[res$ID], function(gs)
    gseaScores(
      geneSet = gs,
      geneList = geneList,
      exponent = exponent
    ))
  if (verbose)
    message("leading edge analysis...")
  ledge <- leading_edge(observed_info)
  res$rank <- ledge$rank
  res$leading_edge <- ledge$leading_edge
  res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse = "/")
  if (verbose)
    message("done...")
  new(
    "gseaResult",
    result = res,
    geneSets = geneSets,
    geneList = geneList,
    params = params,
    readable = FALSE
  )
}


save_gsea <- function(cluster, gsea, path_name, res) {
  gsea_plot <- ggplot(gsea@result, aes(reorder(ID, NES), NES)) +
    geom_col(aes(fill = NES > 0)) +
    scale_fill_manual(values = c("#a6a6a6", "#5a5a5a")) +
    coord_flip() +
    theme_bw(
      base_size = 22,
      base_line_size = 0,
      base_rect_size = 0
    ) +
    theme(legend.position = "none") +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score",
      title = paste0(path_name, " Pathway Enrichment")
    )
  print(gsea_plot)
  ggsave(
    filename = paste0(
      plots_path,
      cluster,
      "_",
      path_name,
      "_chart_res",
      res,
      suffix,
      ".png"
    ),
    height = 14,
    width = 20,
    plot = gsea_plot
  )
  ggsave(
    filename = paste0(
      plots_path,
      cluster,
      "_",
      path_name,
      "_chart_res",
      res,
      suffix,
      ".svg"
    ),
    height = 14,
    width = 20,
    plot = gsea_plot
  )
  
  if (dim(gsea@result)[1] > 1) {
    eplot <- emapplot(gsea, color = "NES", showCategory = 100) +
      scale_color_continuous(low = "blue", high = "red")
    print(eplot)
    ggsave(
      filename = paste0(
        plots_path,
        cluster,
        "_",
        path_name,
        "_eplot_res",
        res,
        suffix,
        ".svg"
      ),
      height = 14,
      width = 20,
      plot = eplot
    )
    ggsave(
      filename = paste0(
        plots_path,
        cluster,
        "_",
        path_name,
        "_eplot_res",
        res,
        suffix,
        ".png"
      ),
      height = 14,
      width = 20,
      plot = eplot
    )
  }
}