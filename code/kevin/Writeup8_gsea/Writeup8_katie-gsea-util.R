calculate_qvalue <- getAnywhere(calculate_qvalue)[["objs"]][[1]]
TERM2NAME <- getAnywhere(TERM2NAME)[["objs"]][[1]]
gseaScores <- getAnywhere(gseaScores)[["objs"]][[1]]
leading_edge <- getAnywhere(leading_edge)[["objs"]][[1]]
build_Anno <- getAnywhere(build_Anno)[["objs"]][[1]]

preparePathwaysAndStats <- getAnywhere(preparePathwaysAndStats)[["objs"]][[1]]

D_fmatch <- getAnywhere(C_fmatch)
fmatch <- function (x, table, nomatch = NA_integer_, incomparables = NULL) {
  .Call(D_fmatch$objs[[1]]$address, x, table, nomatch, incomparables, FALSE)
}

setUpBPPARAM <- getAnywhere(setUpBPPARAM)[["objs"]][[1]]
fgseaSimpleImpl <- getAnywhere(fgseaSimpleImpl)[["objs"]][[1]]
SerialParam <- getAnywhere(SerialParam)[["objs"]][[1]]
multilevelImpl <- function (multilevelPathwaysList, stats, sampleSize, seed, eps, 
                            sign = FALSE, BPPARAM = NULL) {
  size = ES = NULL
  res <- lapply(multilevelPathwaysList, FUN = function(x) {
    fgseaMultilevelCpp(x[, ES], stats, unique(x[, size]),
                       sampleSize, seed, eps, sign)
  })
  return(res)
}
fgseaMultilevelCpp <- getAnywhere(fgseaMultilevelCpp)[["objs"]][[1]]