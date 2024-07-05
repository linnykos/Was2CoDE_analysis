fgsea <- function (pathways,
                   stats,
                   sampleSize = 101,
                   minSize = 1,
                   maxSize = Inf,
                   eps = eps,
                   scoreType = c("std", "pos", "neg"),
                   nproc = 0,
                   gseaParam = 1,
                   BPPARAM = NULL,
                   absEps = NULL)
{
  scoreType <- match.arg(scoreType)
  pp <- preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, scoreType)
  pathwaysFiltered <- pp$filtered
  pathwaysSizes <- pp$sizes
  stats <- pp$stats
  m <- length(pathwaysFiltered)
  if (m == 0) {
    return(
      data.table(
        pathway = character(),
        pval = numeric(),
        padj = numeric(),
        log2err = numeric(),
        ES = numeric(),
        NES = numeric(),
        size = integer(),
        leadingEdge = list()
      )
    )
  }
  if (!is.null(absEps)) {
    warning(
      "You are using deprecated argument `absEps`. ",
      "Use `eps` argument instead. ",
      "`absEps` was assigned to `eps`."
    )
    eps <- absEps
  }
  if (sampleSize < 3) {
    warning("sampleSize is too small, so sampleSize = 3 is set.")
    sampleSize <- max(3, sampleSize)
  }
  log2err = nMoreExtreme = pathway = pval = padj = NULL
  nLeZero = nGeZero = leZeroMean = geZeroMean = nLeEs = nGeEs = isCpGeHalf = NULL
  ES = NES = size = leadingEdge = NULL
  . = "damn notes"
  nPermSimple <- 1000
  minSize <- max(minSize, 1)
  eps <- max(0, min(1, eps))
  if (sampleSize %% 2 == 0) {
    sampleSize <- sampleSize + 1
  }
  gseaStatRes <- do.call(
    rbind,
    lapply(
      pathwaysFiltered,
      calcGseaStat,
      stats = stats,
      returnLeadingEdge = TRUE,
      scoreType = scoreType
    )
  )
  leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
  pathwayScores <- unlist(gseaStatRes[, "res"])
  seeds <- 42
  BPPARAM <- setUpBPPARAM(nproc = nproc, BPPARAM = BPPARAM)
  simpleFgseaRes <- fgseaSimpleImpl(
    pathwayScores = pathwayScores,
    pathwaysSizes = pathwaysSizes,
    pathwaysFiltered = pathwaysFiltered,
    leadingEdges = leadingEdges,
    permPerProc = nPermSimple,
    seeds = seeds,
    toKeepLength = m,
    stats = stats,
    BPPARAM = SerialParam(),
    scoreType = scoreType
  )
  simpleFgseaRes[, `:=`(leZeroMean, NULL)]
  simpleFgseaRes[, `:=`(geZeroMean, NULL)]
  simpleFgseaRes[, `:=`(nLeEs, NULL)]
  simpleFgseaRes[, `:=`(nGeEs, NULL)]
  simpleFgseaRes[, `:=`(nLeZero, NULL)]
  simpleFgseaRes[, `:=`(nGeZero, NULL)]
  unbalanced <- simpleFgseaRes[is.na(pval)]
  unbalanced[, `:=`(padj, as.numeric(NA))]
  unbalanced[, `:=`(log2err, as.numeric(NA))]
  if (nrow(unbalanced) > 0) {
    warning(
      "There were ",
      paste(nrow(unbalanced)),
      " pathways for which P-values were not calculated properly due to ",
      "unbalanced (positive and negative) gene-level statistic values."
    )
  }
  simpleFgseaRes <- simpleFgseaRes[!is.na(pval)]
  simpleError <- 1 / log(2) * sqrt(trigamma(simpleFgseaRes$nMoreExtreme +
                                              1) - trigamma(nPermSimple + 1))
  multError <- sapply((simpleFgseaRes$nMoreExtreme + 1) / (nPermSimple +
                                                             1),
                      multilevelError,
                      sampleSize)
  if (all(multError >= simpleError)) {
    simpleFgseaRes[, `:=`(log2err, 1 / log(2) * sqrt(trigamma(nMoreExtreme +
                                                                1) - trigamma((nPermSimple + 1))))]
    simpleFgseaRes <- rbindlist(list(simpleFgseaRes, unbalanced), use.names = TRUE)
    setorder(simpleFgseaRes, pathway)
    simpleFgseaRes[, `:=`("nMoreExtreme", NULL)]
    setcolorder(
      simpleFgseaRes,
      c(
        "pathway",
        "pval",
        "padj",
        "log2err",
        "ES",
        "NES",
        "size",
        "leadingEdge"
      )
    )
    simpleFgseaRes <- simpleFgseaRes[]
    return(simpleFgseaRes)
  }
  dtSimpleFgsea <- simpleFgseaRes[multError >= simpleError]
  dtSimpleFgsea[, `:=`(log2err, 1 / log(2) * sqrt(trigamma(nMoreExtreme +
                                                             1) - trigamma(nPermSimple + 1)))]
  dtMultilevel <- simpleFgseaRes[multError < simpleError]
  multilevelPathwaysList <- split(dtMultilevel, by = "size")
  indxs <- sample(1:length(multilevelPathwaysList))
  multilevelPathwaysList <- multilevelPathwaysList[indxs]
  
  seed <- 42
  sign <- if (scoreType %in% c("pos", "neg"))
    TRUE
  else
    FALSE
  cpp.res <- multilevelImpl(
    multilevelPathwaysList,
    stats,
    sampleSize,
    seed,
    eps,
    sign = sign,
    BPPARAM = BPPARAM
  )
  cpp.res <- rbindlist(cpp.res)
  result <- rbindlist(multilevelPathwaysList)
  result[, `:=`(pval, cpp.res$cppMPval)]
  result[, `:=`(isCpGeHalf, cpp.res$cppIsCpGeHalf)]
  result[, `:=`(log2err, multilevelError(pval, sampleSize = sampleSize))]
  result[isCpGeHalf == FALSE, `:=`(log2err, NA)]
  if (!all(result$isCpGeHalf)) {
    warning(
      "For some of the pathways the P-values were likely overestimated. ",
      "For such pathways log2err is set to NA."
    )
  }
  result[, `:=`(isCpGeHalf, NULL)]
  result <- rbindlist(list(result, dtSimpleFgsea, unbalanced), use.names = TRUE)
  result[, `:=`(nMoreExtreme, NULL)]
  result[pval < eps, `:=`(c("pval", "log2err"), list(eps, NA))]
  result[, `:=`(padj, p.adjust(pval, method = "BH"))]
  if (nrow(result[pval == eps & is.na(log2err)])) {
    warning(
      "For some pathways, in reality P-values are less than ",
      paste(eps),
      ". You can set the `eps` argument to zero for better estimation."
    )
  }
  setcolorder(result,
              c(
                "pathway",
                "pval",
                "padj",
                "log2err",
                "ES",
                "NES",
                "size",
                "leadingEdge"
              ))
  setorder(result, pathway)
  result <- result[]
  result
}
