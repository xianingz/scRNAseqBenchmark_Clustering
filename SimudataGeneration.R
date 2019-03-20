#Simulation Data Generation

## Define parameter spaces
dropouts <- c(2:6)*0.5
libsizes <- c(10:15)*0.5
libscales <- c(1:5)*0.1

groupequal <- rep(1/4,4)
groupunequal1 <- c(1/3,1/3,1/6,1/6)
groupunequal2 <- c(3/8,1/4,1/4,1/8)
groupmatrix <- cbind(groupequal, groupunequal1, groupunequal2)

distequal <- c(0.005, 0.0075, 0.01, 0.0125, 0.015)
distunequal1 <- c(0.0025,0.0025,0.005,0.01) # * 1:5
distunequal2 <- c(0.0025,0.0025,0.005,0.005) # * 1:5

distmatrix <- matrix(0,nrow=4,ncol=15)
for(i in c(1:5)){
  distmatrix[,i] = rep(distequal[i],4)
}
for(i in (1:5)){
  distmatrix[,i+5] = distunequal1 * i
}
for(i in (1:5)){
  distmatrix[,i+10] = distunequal2 * i
}

#define code
clusterdist <- c(paste0("EDA",c(1:5)), paste0("UDB",c(1:5)),paste0("UDC",c(1:5)))
clustersize <- c("ESA","USB","USC")
libsize <- paste0("LS",c(1:5))
libscale <- paste0("LC",c(1:5))
dropout <- paste0("DO",c(1:5))

library(splatter)
library(scater)

for(i in c(1:15)){
  for(j in c(1:3)){
    for(m in c(1:5)){
      for(n in c(1:5)){
        for(x in c(1:5)){
          params <- newSplatParams()
          params <- setParam(params, "seed",123)
          params <- setParam(params, "batchCells",400)
          params <- setParam(params, "nGenes",10000)
          params <- setParam(params, "de.prob",distmatrix[,i])
          params <- setParam(params, "group.prob", groupmatrix[,j])
          params <- setParam(params, "lib.loc", libsizes[m])
          params <- setParam(params, "lib.scale", libsizes[n])
          params <- setParam(params, "dropout.mid", dropouts[x])
          sim <- splatSimulate(params,method = "groups")
          name <- paste(clusterdist[i], clustersize[j],libsize[m],libscale[n],dropout[x],sep = "_")
          fname = paste0("~/storage/Github/scRNAseqBenchmark_Clustering/Data/Simulations/",name,".rds")
          saveRDS(sim,file = fname)
          
        }
      }
    }
  }
}
