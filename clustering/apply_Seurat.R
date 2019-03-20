## Apply Seurat

suppressPackageStartupMessages({
  library(Seurat)
})

apply_Seurat <- function(sce, params, resolution) {
      (seed <- round(1e6*runif(1)))
      data <- CreateSeuratObject(raw.data = dat, min.cells = params$min.cells,
                                 min.genes = params$min.genes, project = "scRNAseq", 
                                 display.progress = FALSE) 
      data <- NormalizeData(object = data, normalization.method = "LogNormalize", 
                            scale.factor = 1e4, display.progress = FALSE)
      data <- ScaleData(object = data, display.progress = FALSE)
      data <- RunPCA(object = data, pc.genes = rownames(data@data), do.print = FALSE, 
                     pcs.compute = max(params$dims.use), seed.use = seed)
      data <- FindClusters(object = data, reduction.type = "pca", save.SNN = TRUE, 
                           dims.use = params$dims.use, k.param = 30,
                           resolution = resolution, print.output = 0, 
                           random.seed = seed)
      cluster <- data@ident
      cluster
}
