suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
})

postprocessing <- function(sce_filt){
  sce_filt <- calculateQCMetrics(sce_filt)
  sce_filt <- computeSumFactors(sce_filt, sizes = pmin(ncol(sce_filt), seq(20, 120, 20)), min.mean = 0.5)
  sce_filt <- sce_filt[, sizeFactors(sce_filt) > 0]
  sce_filt <- scater::normalize(sce_filt, exprs_values = "counts", return_log = TRUE)
  sce_filt <- scater::normalize(sce_filt, exprs_values = "counts", return_log = FALSE)
  sce_filt <- runTSNE(sce_filt, exprs_values = "logcounts", perplexity = 10)
  sce_filt
}

