## Apply GiniClust2

source("/home/xianingz/storage/Scbenchmarking/scRNAseqClusteringBenchmarking/Rscripts/clustering/GiniClust2/GiniClust2.R")

apply_GiniClust2 <- function(sce, params, k){
    dat <- counts(sce)
    cluster <- GiniClust2(sce = sce, id = "Simu", k = k, MinPts = params$MinPts, eps = params$eps, gap_statistic = FALSE, automatic_eps = FALSE, automatic_minpts = FALSE)
    names(cluster) <- colnames(dat)
    cluster
    #Determine number of clusters automatically
    #est_k <- length(unique(GiniClust2(sce = sce, id = "Simu")))
}
