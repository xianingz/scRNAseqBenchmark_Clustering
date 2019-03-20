## Apply SIMLR

suppressPackageStartupMessages({
  library(SIMLR)
})

apply_SIMLR <- function(sce, params, k){
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      cluster <- SIMLR(X = dat, c = k)$y$cluster
      names(cluster) <- colnames(dat)
    })
    #Determine number of clusters automatically
    est_k <- SIMLR_Estimate_Number_of_Clusters(X = dat, NUMC = 2:6)
    est_k <- round(mean(c(c(2:6)[which.min(est_k$K1)],c(2:6)[which.min(est_k$K2)])))
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e){
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA),
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
