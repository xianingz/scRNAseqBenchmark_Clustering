##Generate PBS files
library(stringr)
names <- list.files("~/storage/Github/scRNAseqBenchmark_Clustering/Data/data_pre_do/")
#names <- names[grepl("EDA",names)]
names <- str_replace_all(names,".rds","")
cont1 = '#!/bin/sh
####  PBS preamble\n'

cont3 = '# Change "bjensen" to your uniqname:
#PBS -M xianingz@umich.edu
#PBS -m abe

# Change the number of cores (ppn=1), amount of memory, and walltime:
#PBS -l nodes=1:ppn=1,pmem=8gb,walltime=8:00:00
#PBS -j oe
#PBS -V

# Change "example_flux" to the name of your Flux allocation:
#PBS -A junzli_flux
#PBS -q flux
#PBS -l qos=flux

####  End PBS preamble

#  Show list of CPUs you ran on, if you\'re running under PBS
if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi

#  Change to the directory you submitted from
if [ -n "$PBS_O_WORKDIR" ]; then cd $PBS_O_WORKDIR; fi

#  Put your job commands here:
Rscript /home/xianingz/storage/Github/scRNAseqBenchmark_Clustering/benchmark.R '
for(name in names){
  cont2 =paste0('#PBS -N ', name,'\n')
  writeLines(paste0(cont1,cont2, cont3, name), paste0('/home/xianingz/storage/Github/scRNAseqBenchmark_Clustering/pbs_scripts/flux_pbs/',name,".pbs"))
}

