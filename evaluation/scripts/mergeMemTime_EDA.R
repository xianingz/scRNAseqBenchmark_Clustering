library(stringr)
library(SummarizedExperiment)
names <- list.files("~/storage/Github/scRNAseqBenchmark_Clustering/results/pcareduce/")
names <- names[grepl("EDA",names)]
names <- str_replace_all(names,".res.rds","")
mem <- matrix(0, nrow=225, ncol=length(names))
time <- matrix(0, nrow=225, ncol=length(names))
colnames(mem) = names
colnames(time) = names

for(name in names){
  if(grepl("DO4", name)){
    pr <- readRDS(paste0("~/storage/Github/scRNAseqBenchmark_Clustering/results/pcareduce/",name,".res.rds"))
    or <- readRDS(paste0("~/storage/Github/scRNAseqBenchmark_Clustering/results/wopcareduce/",name,".res.rds"))
    newres = or@colData
    newres[rownames(pr@colData),] = pr@colData
    newres <- newres[!grepl("SAFE",rownames(newres)),]
    if(is.null(rownames(mem))){
      rownames(mem) = rownames(newres)
      rownames(time) = rownames(newres)
      mem[,name] = newres$mem_peak
      time[,name] = newres$time_change
    }else{
      mem[rownames(newres),name] = newres$mem_peak
      time[rownames(newres),name] = newres$time_change
    }
  }else{
    rds <- readRDS(paste0("~/storage/Github/scRNAseqBenchmark_Clustering/results/", name, ".res.rds"))
    newres = rds@colData
    if(is.null(rownames(mem))){
      rownames(mem) = rownames(newres)
      rownames(time) = rownames(newres)
      mem[,name] = newres$mem_peak
      time[,name] = newres$time_change
    }else{
      mem[rownames(newres),name] = newres$mem_peak
      time[rownames(newres),name] = newres$time_change
    }
  }
}
saveRDS(mem, file = paste0("~/storage/Github/scRNAseqBenchmark_Clustering/evaluation/mem.EDA.rds"))
saveRDS(time, file = paste0("~/storage/Github/scRNAseqBenchmark_Clustering/evaluation/time.EDA.rds"))
