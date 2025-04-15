#' UMAP Maker
#'
#' This function generates UMAP embeddings from PCA data.
#'
#' @param PCA_data A data frame containing PCA results with a column 'u' containing PCA coordinates.
#' @param n_neighbors Integer specifying the number of neighbors to consider for UMAP.
#' @param metric Character string indicating the distance metric to use (default is "correlation").
#'
#' @return A matrix of UMAP results.
#' 
#' @export
#' @import uwot
Umap_maker = function(PCA_data,n_neighbors =30, metric = "correlation"){
  Umap_result <- uwot::umap(PCA_data$u, n_neighbors = n_neighbors, metric = metric, verbose = TRUE)
  return(Umap_result)
}

