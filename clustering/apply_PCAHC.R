## Apply PCA + Hierarchical clustering

apply_PCAHC <- function(sce, params, k) {
  dat <- logcounts(sce)
  pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
  pca <- pca$x[, seq_len(params$nPC), drop = FALSE]
  hcl <- hclust(dist(pca), method = "ward.D2")
  cluster <- cutree(hcl, k = k)
  names(cluster) = colnames(dat)
  cluster
}