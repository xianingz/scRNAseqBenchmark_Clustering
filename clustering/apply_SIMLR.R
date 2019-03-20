## Apply SIMLR

suppressPackageStartupMessages({
  library(SIMLR)
})

apply_SIMLR <- function(sce, params, k){
    dat <- logcounts(sce)
    cluster <- SIMLR(X = dat, c = k)$y$cluster
    names(cluster) <- colnames(dat)
    cluster
    #Determine number of clusters automatically
    #est_k <- SIMLR_Estimate_Number_of_Clusters(X = dat, NUMC = 2:6)
    #est_k <- round(mean(c(c(2:6)[which.min(est_k$K1)],c(2:6)[which.min(est_k$K2)])))
}
