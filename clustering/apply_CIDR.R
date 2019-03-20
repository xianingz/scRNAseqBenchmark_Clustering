## Apply CIDR

suppressPackageStartupMessages({
  library(cidr)
})

apply_CIDR <- function(sce, params, k) {
      dat <- counts(sce)
      sData <- scDataConstructor(dat, tagType = "raw")
      sData <- determineDropoutCandidates(sData)
      sData <- wThreshold(sData)
      sData <- scDissim(sData, threads = 1)
      sData <- scPCA(sData, plotPC = FALSE)
      sData <- nPC(sData)
      
      ## Cluster with preset number of clusters
      sDataC <- scCluster(object = sData, nCluster = k, 
                          nPC = sData@nPC, cMethod = "ward.D2")
      cluster <- sDataC@clusters
      names(cluster) <- colnames(sDataC@tags)
      cluster
}
