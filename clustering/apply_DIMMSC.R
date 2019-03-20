## Apply SIMLR

suppressPackageStartupMessages({
  library(DIMMSC)
})

apply_DIMMSC <- function(sce, params, k){
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      cluster <- DIMMSC(data=dat, K=k, method_cluster_intial="kmeans", method_alpha_intial="Ronning", maxiter=200, tol=1e-4, lik.tol=1e-2)$mem
      names(cluster) <- colnames(dat)
    })
    #Cannot determine number of clusters automatically
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e){
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA),
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
