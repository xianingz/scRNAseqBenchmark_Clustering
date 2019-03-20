## Apply FlowSOM

suppressPackageStartupMessages({
  library(flowCore)
  library(FlowSOM)
})

apply_FlowSOM <- function(sce, params, k) {
  dat <- logcounts(sce)
  pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
  pca <- pca$x[, seq_len(params$nPC), drop = FALSE]
  ff <- flowFrame(exprs = pca)
  fSOM <- FlowSOM::ReadInput(ff, compensate = FALSE, transform = FALSE, 
                                 scale = FALSE, silent = TRUE)
  fSOM <- FlowSOM::BuildSOM(fSOM, silent = TRUE, xdim = params$xdim, 
                                ydim = params$ydim)
  metaClustering <- metaClustering_consensus(fSOM$map$codes, k = k)
  cluster <- metaClustering[fSOM$map$mapping[, 1]]
  names(cluster) <- colnames(dat)
  cluster
}
