filterExpr <- function(sce, pctkeep) {
  exprsn <- rowMeans(logcounts(sce))
  keep <- order(exprsn, decreasing = TRUE)[seq_len(pctkeep/100*length(exprsn))]
  return(postprocessing(sce[keep, ]))
}
