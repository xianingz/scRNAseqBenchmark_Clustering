GiniClust2 <- function(sce, id, gap_statistic=TRUE, automatic_eps=TRUE, automatic_minpts=TRUE, k = NULL, MinPts = NULL, eps = NULL){
  minCellNum = 0 #default 3
  minGeneNum = 0 #default 2000
  expressed_cutoff = 1
  gini.bi = 0
  log2.expr.cutoffl = 0
  log2.expr.cutoffh = 20
  Gini.pvalue_cutoff = 0.0001
  Norm.Gini.cutoff = 1
  span = 0.9
  outlier_remove = 0.75
  GeneList = 1
  Gamma = 0.9
  diff.cutoff = 1
  lr.p_value_cutoff = 1e-5
  CountsForNormalized = 100000
  
  #Paths
  workdir = "/home/xianingz/storage/Scbenchmarking/scRNAseqClusteringBenchmarking/Rscripts/clustering/GiniClust2/GiniClust2_download/Proj/Simulation/"
  Rfundir = "/home/xianingz/storage/Scbenchmarking/scRNAseqClusteringBenchmarking/Rscripts/clustering/GiniClust2/GiniClust2_download/Rfunction/"
  exprimentID = "simu"
  
  setwd(workdir)
  dir.create(file.path(workdir, "results"), showWarnings = FALSE)
  dir.create(file.path(workdir, "figures"), showWarnings = FALSE)
  
  #Dataset-specific parameters
  #MinPts = 3
  #eps = 0.5
  #k = 3
  #gap_statistic = TRUE
  K.max = 10
  #automatic_eps = TRUE
  #automatic_minpts = TRUE
  perplexity_G = 30
  perplexity_F = 30
  max_iter_G = 1000
  max_iter_F = 1000
  
  #Load packages and functions
  source(paste(Rfundir, "GiniClust2_packages.R",sep = ""), local=TRUE)
  source(paste(Rfundir, "GiniClust2_functions.R",sep = ""), local=TRUE)
  
  #Load, preprocess and filter data
  data <- counts(sce)
  source(paste(Rfundir,"GiniClust2_preprocess.R",sep = ""), local=TRUE)
  source(paste(Rfundir,"GiniClust2_filtering_RawCounts.R",sep = ""), local=TRUE)
  
  #Run GiniClust
  library(dbscan)
  source(paste(Rfundir,"GiniClust2_fitting.R",sep = ""), local=TRUE)
  source(paste(Rfundir,"GiniClust2_Gini_clustering.R",sep = ""), local=TRUE)
  #table(P_G)
  #source(paste(Rfundir,"GiniClust2_Gini_tSNE.R",sep = ""))
  
  #Run Fano-based k-means clustering steps
  source(paste(Rfundir,"GiniClust2_Fano_clustering.R",sep = ""), local=TRUE)
  #table(P_F)
  #source(paste(Rfundir,"GiniClust2_Fano_tSNE.R",sep = ""))
  
  #Run cluster-aware weighted consensus clustering
  source(paste(Rfundir,"GiniClust2_consensus_clustering.R",sep = ""), local=TRUE)
  return(finalCluster)
  
}
