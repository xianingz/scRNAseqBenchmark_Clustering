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

#Generate the code for expand methods
methods <- c("SC3","CIDR","DIMMSC","FlowSOM","GiniClust2","monocle","PCAHC",
             "PCAKmeans","pcaReduce","RaceID2","RtsneKmeans","SAFE","SC3svm","SIMLR","TSCAN")
for(i in methods){
  cat('  expandMethodInWorkflow(step = "Clustering", label = "',i,'", 
                         onlyone = "k",
                         params = rlang::quos(',i,'_2 = 2,\n',
      '                                           ',i,'_3 = 3,\n',
      '                                           ',i,'_5 = 5,\n',
      '                                           ',i,'_6 = 6)) %>%\n', sep="")
}