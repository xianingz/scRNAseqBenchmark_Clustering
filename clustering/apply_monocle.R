## Apply monocle

suppressPackageStartupMessages({
  library(monocle)
  library(scran)
})

apply_monocle <- function(sce, params, k) {
    cds <- convertTo(sce, type = "monocle")
    cds <- tryCatch({
      estimateDispersions(cds)
    }, error = function(e) {
      cds
    })
    cds <- reduceDimension(cds, max_components = params$max_components, 
                             num_dim = params$num_dim,
                             reduction_method = "tSNE", verbose = TRUE)
    cds <- clusterCells(cds, num_clusters = k + 1, method = "densityPeak")
    cluster <- cds$Cluster
    names(cluster) <- colnames(cds)
    cluster
}
