## Apply TSCAN

suppressPackageStartupMessages({
  library(TSCAN)
})

apply_TSCAN <- function(sce, params, k) {
    dat <- logcounts(sce)
    ## Remove genes with variance = 0
    dat <- dat[rowVars(dat) > 0, ]
    cluster <- exprmclust(dat, clusternum = k, modelNames = "VVV", reduce = TRUE)$clusterid
    cluster
    ## Determine number of clusters automatically
    #est_k <- length(unique(exprmclust(dat, clusternum = params$range_clusters, 
    #                                  modelNames = "VVV", reduce = TRUE)$clusterid))
    
}
