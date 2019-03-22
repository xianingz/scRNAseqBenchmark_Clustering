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
source("clustering/apply_SAFE.R")

##test filter function
sce <- readRDS("Data/data_pre/EDA1_ESA_LS3_LC4_DO4.rds")
sce.f <- filterHVG(sce, pctkeep = 20)
sce.f <- filterExpr(sce, pctkeep = 20)
sce.f <- filterM3Drop(sce, pctkeep = 20)

##test cluster function
params <- list()
cluster <- apply_CIDR(sce.f, params = params, k=4)

params <- list()
cluster <- apply_DIMMSC(sce.f, params = params, k=4)

params <- list(nPC=30, xdim=5, ydim=5)
cluster <- apply_FlowSOM(sce.f,params = params, k=4)

params <- list(MinPts=3, eps=0.5)
cluster <- apply_GiniClust2(sce, params = params,k=4)

params <- list(num_dim = 50, max_components = 3)
cluster <- apply_monocle(sce.f, params = params, k=4)

params <- list(nPC=30)
cluster <- apply_PCAHC(sce.f, params = params, k=4)

params <- list(nPC=30)
cluster <- apply_PCAKmeans(sce.f, params = params, k=4)

params <- list(nbt = 100, q = 30)
cluster <- apply_pcaReduce(sce.f, params = params, k=4)

params <- list(mintotal = 1, minexpr = 0, minnumber = 1, maxexpr = Inf, range_clusters=2:6)
cluster <- apply_RaceID2(sce.f, params = params, k=4)

params <- list(perplexity = 30, dims = 3, initial_dims = 50)
cluster <- apply_RtsneKmeans(sce.f, params = params, k=4) 

params <- list()
cluster <- apply_SAFE(sce.f, params = params, k=4)

params <- list(pct_droupout_min = 0, pca_dropout_max = 100, gene_filter = FALSE)
cluster <- apply_SC3(sce.f, params = params, k=4)

params <- list(pct_droupout_min = 0, pca_dropout_max = 100, gene_filter = FALSE)
cluster <- apply_SC3svm(sce.f, params = params, k=4)

params <- list(min.cells = 0, min.genes = 0, dims.use = 1:30)
cluster <- apply_Seurat(sce.f, params = params, resolution = 1)

params <- list()
cluster <- apply_SIMLR(sce.f, params = params, k=4)

params <- list()
cluster <- apply_TSCAN(sce.f, params = params, k=4)

library(SummarizedBenchmark)
source("")
