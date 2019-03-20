## Apply GiniClust2

source("/home/xianingz/storage/Scbenchmarking/scRNAseqClusteringBenchmarking/Rscripts/clustering/GiniClust2/GiniClust2.R")

apply_GiniClust2 <- function(sce, params, k){
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      cluster <- GiniClust2(sce = sce, id = "Simu", k = k, MinPts = params$MinPts, eps = params$eps, gap_statistic = FALSE, automatic_eps = FALSE, automatic_minpts = FALSE)
      names(cluster) <- colnames(dat)
    })
    #Determine number of clusters automatically
    est_k <- length(unique(GiniClust2(sce = sce, id = "Simu")))
    
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
