## Apply SIMLR

suppressPackageStartupMessages({
  library(DIMMSC)
})

apply_DIMMSC <- function(sce, params, k){
    dat <- counts(sce)
    cluster <- DIMMSC(data=dat, K=k, method_cluster_intial="kmeans", method_alpha_intial="Ronning", maxiter=200, tol=1e-4, lik.tol=1e-2)$mem
    names(cluster) <- colnames(dat)
    cluster
}
