library(stringr)
library(SummarizedExperiment)
names <- list.files("~/storage/Github/scRNAseqBenchmark_Clustering/results/pcareduce/")
names <- names[grepl("EDA2",names)]
names <- str_replace_all(names,".res.rds","")
res <- list()
tru <- list()

listtomatrix <- function(li){
  unli <- apply(li,2,unlist)
  ma <- matrix(0,nrow = length(unli[[1]]), ncol=dim(li)[2])
  for(i in c(1:dim(li)[2])){
    if(length(unli[[i]]) == 1){
      ma[,i] = rep(NA, length(unli[[1]]))
    }else if(length(unli[[i]]) > length(unli[[1]])){
      ma[,i] = unli[[i]][1:length(unli[[1]])]
    }else{
      ma[,i] = unli[[i]]
    }
  }
  ma
}

for(name in names){
  if(grepl("DO4", name)){
    pr <- readRDS(paste0("~/storage/Github/scRNAseqBenchmark_Clustering/results/pcareduce/",name,".res.rds"))
    or <- readRDS(paste0("~/storage/Github/scRNAseqBenchmark_Clustering/results/wopcareduce/",name,".res.rds"))
    newres = assay(or)
    if(dim(assay(or))[1] == 1){
      newres = listtomatrix(newres)
    }
    colnames(newres) = colnames(assay(or))
    pcres = assay(pr)
    if(dim(assay(pr))[1] == 1){
      pcres = listtomatrix(pcres)
    }  
    newres[,colnames(assay(pr))] = pcres
    newres <- newres[,!grepl("SAFE",colnames(newres))]
    res[[name]] = newres
    tru[[name]] = metadata(pr)$sessions[[1]]$parameters$truthCols
  }else{
    rds <- readRDS(paste0("~/storage/Github/scRNAseqBenchmark_Clustering/results/", name, ".res.rds"))
    rdsmat <- assay(rds)
    if(dim(rdsmat)[1]==1){
      rdsmat = listtomatrix(rdsmat)
    }
    res[[name]] <- rdsmat
    tru[[name]] <- metadata(rds)$sessions[[1]]$parameters$truthCols
  }
}
saveRDS(res, file = paste0("~/storage/Github/scRNAseqBenchmark_Clustering/evaluation/res.EDA2.rds"))
saveRDS(tru, file = paste0("~/storage/Github/scRNAseqBenchmark_Clustering/evaluation/tru.EDA2.rds"))
