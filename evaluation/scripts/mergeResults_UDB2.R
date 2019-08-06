##merge all datasets
library(stringr)
library(SummarizedExperiment)
names <- list.files("~/storage/Github/scRNAseqBenchmark_Clustering/results/")
names <- names[grepl("UDB2",names)]
#names <- str_replace_all(names,".res.rds","")
res <- list()
tru <- list()
i=1

listtomatrix <- function(li){
  unli <- apply(li,2,unlist)
  ma <- matrix(0,nrow = length(unli[[1]]), ncol=225)
  for(i in c(1:225)){
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
  cat(i)
  i = i + 1
  rds <- readRDS(paste0("~/storage/Github/scRNAseqBenchmark_Clustering/results/", name))
  name <- str_replace(name,".res.rds","")
  rdsmat <- assay(rds)
  if(dim(rdsmat)[1]==1){
    rdsmat = listtomatrix(rdsmat)
  }
  res[[name]] <- rdsmat
  tru[[name]] <- metadata(rds)$sessions[[1]]$parameters$truthCols
}

saveRDS(res, file = "~/storage/Github/scRNAseqBenchmark_Clustering/evaluation/res.UDB2.rds")
saveRDS(tru, file = "~/storage/Github/scRNAseqBenchmark_Clustering/evaluation/tru.UDB2.rds")
