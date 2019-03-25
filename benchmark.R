#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(R.utils)
library(SummarizedBenchmark)
sourceDirectory("~/storage/Github/SummarizedBenchmark_Workflow/R/",modifiedOnly=FALSE)

###Read in gene filter functions
setwd("~/storage/Github/scRNAseqBenchmark_Clustering/")
source("geneselection/filter_genes.R")
source("geneselection/filterExpr.R")
source("geneselection/filterHVG.R")
source("geneselection/filterM3Drop.R")

###Read in cluster functions
source("clustering/apply_CIDR.R")
source("clustering/apply_DIMMSC.R")
source("clustering/apply_FlowSOM.R")
source("clustering/apply_GiniClust2.R")
source("clustering/apply_monocle.R")
source("clustering/apply_PCAHC.R")
source("clustering/apply_PCAKmeans.R")
source("clustering/apply_pcaReduce.R")
source("clustering/apply_RaceID2.R")
source("clustering/apply_RtsneKmeans.R")
source("clustering/apply_TSCAN.R")
source("clustering/apply_SIMLR.R")
source("clustering/apply_Seurat.R")
source("clustering/apply_SC3svm.R")
source("clustering/apply_SC3.R")

##test filter function
sce <- readRDS(paste0("/home/xianingz/storage/Github/scRNAseqBenchmark_Clustering/Data/data_pre/",args[1],".rds"))
b <- BenchDesign(data = list(tdata = sce))
b <- initWorkflow(b, steps=c("Genefilter","Clustering"))
b <- b %>% 
  addMethodToWorkflow(step="Genefilter",
                      method_label = "HVG",
                      func = filterHVG,
                      params = rlang::quos(sce = tdata, pctkeep = 20)) %>%
  addMethodToWorkflow(step="Genefilter",
                      method_label = "HEG",
                      func = filterExpr,
                      params = rlang::quos(sce = tdata, pctkeep = 20)) %>%
  addMethodToWorkflow(step="Genefilter",
                      method_label = "HDG",
                      func = filterM3Drop,
                      params = rlang::quos(sce = tdata, pctkeep = 20)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "Seurat",
                      func = apply_Seurat,
                      params = rlang::quos(params = list(min.cells = 0, min.genes = 0, dims.use = 1:30), 
                                           resolution = 1)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "SC3",
                      func = apply_SC3,
                      params = rlang::quos(params = list(pct_droupout_min = 0, pca_dropout_max = 100, gene_filter = FALSE),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "CIDR",
                      func = apply_CIDR,
                      params = rlang::quos(params = list(),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "DIMMSC",
                      func = apply_DIMMSC,
                      params = rlang::quos(params = list(),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "FlowSOM",
                      func = apply_FlowSOM,
                      params = rlang::quos(params = list(nPC=30, xdim=5, ydim=5),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "GiniClust2",
                      func = apply_GiniClust2,
                      params = rlang::quos(params = list(MinPts=3, eps=0.5),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "monocle",
                      func = apply_monocle,
                      params = rlang::quos(params = list(num_dim = 50, max_components = 3),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "PCAHC",
                      func = apply_PCAHC,
                      params = rlang::quos(params = list(nPC=30),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "PCAKmeans",
                      func = apply_PCAKmeans,
                      params = rlang::quos(params = list(nPC=30),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "pcaReduce",
                      func = apply_pcaReduce,
                      params = rlang::quos(params = list(nbt = 100, q = 30),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "RaceID2",
                      func = apply_RaceID2,
                      params = rlang::quos(params = list(mintotal = 1, minexpr = 0, minnumber = 1, maxexpr = Inf, range_clusters=2:6),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "RtsneKmeans",
                      func = apply_RtsneKmeans,
                      params = rlang::quos(params = list(perplexity = 30, dims = 3, initial_dims = 50),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "SC3svm",
                      func = apply_SC3svm,
                      params = rlang::quos(params = list(pct_droupout_min = 0, pca_dropout_max = 100, gene_filter = FALSE),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "SIMLR",
                      func = apply_SIMLR,
                      params = rlang::quos(params = list(),
                                           k = 4)) %>%
  addMethodToWorkflow(step="Clustering",
                      method_label = "TSCAN",
                      func = apply_TSCAN,
                      params = rlang::quos(params = list(),
                                           k = 4))

## expand parameter space
b <- b %>% 
  expandMethodInWorkflow(step = "Clustering", label = "Seurat", 
                         onlyone = "resolution",
                         params = rlang::quos(Seurat_2 = 0.6,
                                              Seurat_3 = 0.8,
                                              Seurat_5 = 1.2,
                                              Seurat_6 = 1.4)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "SC3", 
                         onlyone = "k",
                         params = rlang::quos(SC3_2 = 2,
                                              SC3_3 = 3,
                                              SC3_5 = 5,
                                              SC3_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "CIDR", 
                         onlyone = "k",
                         params = rlang::quos(CIDR_2 = 2,
                                              CIDR_3 = 3,
                                              CIDR_5 = 5,
                                              CIDR_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "DIMMSC", 
                         onlyone = "k",
                         params = rlang::quos(DIMMSC_2 = 2,
                                              DIMMSC_3 = 3,
                                              DIMMSC_5 = 5,
                                              DIMMSC_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "FlowSOM", 
                         onlyone = "k",
                         params = rlang::quos(FlowSOM_2 = 2,
                                              FlowSOM_3 = 3,
                                              FlowSOM_5 = 5,
                                              FlowSOM_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "GiniClust2", 
                         onlyone = "k",
                         params = rlang::quos(GiniClust2_2 = 2,
                                              GiniClust2_3 = 3,
                                              GiniClust2_5 = 5,
                                              GiniClust2_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "monocle", 
                         onlyone = "k",
                         params = rlang::quos(monocle_2 = 2,
                                              monocle_3 = 3,
                                              monocle_5 = 5,
                                              monocle_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "PCAHC", 
                         onlyone = "k",
                         params = rlang::quos(PCAHC_2 = 2,
                                              PCAHC_3 = 3,
                                              PCAHC_5 = 5,
                                              PCAHC_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "PCAKmeans", 
                         onlyone = "k",
                         params = rlang::quos(PCAKmeans_2 = 2,
                                              PCAKmeans_3 = 3,
                                              PCAKmeans_5 = 5,
                                              PCAKmeans_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "pcaReduce", 
                         onlyone = "k",
                         params = rlang::quos(pcaReduce_2 = 2,
                                              pcaReduce_3 = 3,
                                              pcaReduce_5 = 5,
                                              pcaReduce_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "RaceID2", 
                         onlyone = "k",
                         params = rlang::quos(RaceID2_2 = 2,
                                              RaceID2_3 = 3,
                                              RaceID2_5 = 5,
                                              RaceID2_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "RtsneKmeans", 
                         onlyone = "k",
                         params = rlang::quos(RtsneKmeans_2 = 2,
                                              RtsneKmeans_3 = 3,
                                              RtsneKmeans_5 = 5,
                                              RtsneKmeans_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "SC3svm", 
                         onlyone = "k",
                         params = rlang::quos(SC3svm_2 = 2,
                                              SC3svm_3 = 3,
                                              SC3svm_5 = 5,
                                              SC3svm_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "SIMLR", 
                         onlyone = "k",
                         params = rlang::quos(SIMLR_2 = 2,
                                              SIMLR_3 = 3,
                                              SIMLR_5 = 5,
                                              SIMLR_6 = 6)) %>%
  expandMethodInWorkflow(step = "Clustering", label = "TSCAN", 
                         onlyone = "k",
                         params = rlang::quos(TSCAN_2 = 2,
                                              TSCAN_3 = 3,
                                              TSCAN_5 = 5,
                                              TSCAN_6 = 6))


b <- buildMethodFromWorkflow(b)
sb <- buildBench(b, truthCols = sce$Group)
saveRDS(sb, file = paste0("/home/xianingz/storage/Github/scRNAseqBenchmark_Clustering/results/",args[1],".res.rds"))
