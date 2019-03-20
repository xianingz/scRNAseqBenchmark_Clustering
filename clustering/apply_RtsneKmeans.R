## Apply t-SNE + K-means. 

suppressPackageStartupMessages({
  library(Rtsne)
})

apply_RtsneKmeans <- function(sce, params, k) {
  dat <- logcounts(sce)
  rtsne <- Rtsne(X = t(dat), dims = params$dims, 
                     perplexity = params$perplexity, pca = TRUE, 
                     initial_dims = params$initial_dims, check_duplicates = FALSE)
  cluster <- kmeans(rtsne$Y, centers = k, nstart = 25)$cluster
  names(cluster) = colnames(dat)
  cluster
}
