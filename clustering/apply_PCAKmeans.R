## Apply PCA + K-means. 

apply_PCAKmeans <- function(sce, params, k) {
  dat <- logcounts(sce)
  pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
  pca <- pca$x[, seq_len(params$nPC), drop = FALSE]
  cluster <- kmeans(pca, centers = k, nstart = 25)$cluster
  names(cluster) = colnames(dat)
  cluster
}