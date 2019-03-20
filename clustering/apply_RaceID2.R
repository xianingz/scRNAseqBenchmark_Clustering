## Apply RaceID2

source("Rscripts/clustering/RaceID2_StemID_class.R")

apply_RaceID2 <- function(sce, params, k) {
  (seed <- round(1e6*runif(1)))
  dat <- as.data.frame(counts(sce))
  sc <- SCseq(dat)
  sc <- filterdata(sc, mintotal = params$mintotal, minexpr = params$minexpr, 
                       minnumber = params$minnumber, maxexpr = params$maxexpr, 
                       downsample = FALSE, dsn = 1, rseed = seed)
      
      ## Cluster with predetermined number of clusters
  cluster <- clustexp(sc, metric = "pearson", cln = k, 
                          do.gap = FALSE, clustnr = max(params$range_clusters), B.gap = 50,
                          SE.method = "Tibs2001SEmax", SE.factor = 0.25, 
                          bootnr = 50, rseed = seed)@cluster$kpart
  cluster
    ## Determine number of clusters automatically
  #est_k <- length(unique(clustexp(sc, metric = "pearson", cln = 0, 
  #                                  do.gap = TRUE, clustnr = max(params$range_clusters), B.gap = 50,
  #                                  SE.method = "Tibs2001SEmax", SE.factor = 0.25, 
  #                                  bootnr = 50, rseed = seed)@cluster$kpart))

}
