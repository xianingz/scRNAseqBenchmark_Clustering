## Apply SAFE clustering

source("Rscripts/clustering/SAFE_2.1_Linux/individual_clustering_modified.R")
source("Rscripts/clustering/SAFE_2.1_Linux/SAFE_modified.R")

apply_SAFE <- function(sce, params, k) {
    dat <- counts(sce)
    indc <- individual_clustering(inputTags = dat, datatype = "count", 
                                    SC3 = TRUE, gene_filter = FALSE, 
                                    svm_num_cells = 1000, CIDR = TRUE, 
                                    Seurat = TRUE, nPC = NULL, 
                                    resolution = 0.9, seurat_min_cell = 0, 
                                    resolution_min = 1.2, tSNE = TRUE, 
                                    var_genes = NULL, SEED = NULL,
                                    n_cores = 1)
    safe <- SAFE(cluster_results = indc, k_min = k, k_max = k, 
                   MCLA = TRUE, HGPA = FALSE, CSPA = FALSE, 
                   cspc_cell_max = NULL)
      
    cluster <- safe$optimal_clustering
    names(cluster) <- colnames(dat)
    cluster
    ## Determine number of clusters automatically
    #safe <- SAFE(cluster_results = indc, k_min = 2, k_max = NULL, 
    #             MCLA = TRUE, HGPA = FALSE, CSPA = FALSE, 
    #             cspc_cell_max = NULL)
    #est_k <- safe$optimal_k
}
