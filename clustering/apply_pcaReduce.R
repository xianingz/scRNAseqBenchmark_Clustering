## Apply pcaReduce.

suppressPackageStartupMessages({
  library(pcaReduce)
  library(clue)
})

apply_pcaReduce <- function(sce, params, k) {
  dat <- logcounts(sce)
  pca <- PCAreduce(t(dat), nbt = params$nbt, q = params$q, method = "S")
  part <- lapply(pca, function(x) {
        colnames(x) <- paste0("k", (params$q + 1):2)
        as.cl_partition(x[, paste0("k", k)])
      })
  cons <- cl_consensus(as.cl_ensemble(part), method = "SE", control = list(nruns = 50))
  cluster <- c(cl_class_ids(cons))
  cluster
}
