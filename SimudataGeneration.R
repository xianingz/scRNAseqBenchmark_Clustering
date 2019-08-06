#Simulation Data Generation

## Define parameter spaces
dropouts <- seq(6,0,-2)*0.5
libsizes <- c(11:14)
libscales <- c(4:1)*0.1

groupequal <- rep(1/4,4)
groupunequal1 <- c(1/3,1/3,1/6,1/6)
groupunequal2 <- c(3/8,1/4,1/4,1/8)
groupmatrix <- cbind(groupequal, groupunequal1, groupunequal2)

#distequal <- c(0.0075, 0.01, 0.0125, 0.015)
distequal <- c(0.01, 0.015, 0.02, 0.025)
distunequal1 <- c(0.01,0.01,0.015,0.015) # * (1,2.5,0.5)
distunequal2 <- c(0.015,0.015,0.01,0.01) # * (1,2.5,0.5)

distmatrix <- matrix(0,nrow=4,ncol=12)
for(i in c(1:4)){
  distmatrix[,i] = rep(distequal[i],4)
}
for(i in (1:4)){
  distmatrix[,i+4] = distunequal1 * (i/4+3/4)
}
for(i in (1:4)){
  distmatrix[,i+8] = distunequal2 * (i/4+3/4)
}

#define code
clusterdist <- c(paste0("EDA",c(1:4)), paste0("UDB",c(1:4)),paste0("UDC",c(1:4)))
clustersize <- c("ESA","USB","USC")
libsize <- paste0("LS",c(1:4))
libscale <- paste0("LC",c(1:4))
dropout <- paste0("DO",c(1:4))

library(splatter)
library(scater)


##3/22/2019 Only for equal cluster distance 4*2*4*4*4 = 512
#for(i in c(1:4)){ ## cluster distance, 4 degrees of difficulty
##3/23/2019 for unequal distance 4*2*4*4*4 = 512
#for(i in c(5:8)){ ## cluster distance, unequal, 4 degrees of difficulty
##6/29/2019 resimulate cluster with dropout.
for(i in c(1:8)){
  for(j in c(1:2)){  ## cluster size, only one unequal is included beacuase it symmetric
    for(m in c(1:4)){  ## library size
      for(n in c(1:4)){  ## library scale
        for(x in c(1:3)){  ## dropout rate
          params <- newSplatParams()
          params <- setParam(params, "seed",123)
          params <- setParam(params, "batchCells",400)
          params <- setParam(params, "nGenes",10000)
          params <- setParam(params, "group.prob", as.vector(groupmatrix[,j]))
          params <- setParam(params, "de.prob", c(distmatrix[,i]))
          params <- setParam(params, "lib.loc", libsizes[m])
          params <- setParam(params, "lib.scale", libscales[n])
          params <- setParam(params, "dropout.type", "experiment")
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

