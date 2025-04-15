#' Perform cell clustering using a K-nearest neighbors approach and Louvain algorithm.
#'
#' @param Data_correction A matrix of reduced expression data.
#' @param K An integer specifying the number of neighbors to consider.
#' @param metric_used A string indicating the distance metric to use.
#' @param nThreads An integer specifying the number of threads for computation.
#' @param resolution A numeric value for the resolution parameter in clustering.
#' @return A character vector indicating the cluster membership for each cell.
#' @export
#' 
Cell_clustering_function = function(Data_correction, K, metric_used, nThreads, resolution) {
  KNN_matrix = Knn(as.matrix(Data_correction), k = K, verbose = TRUE, indexType = metric_used, nThreads = nThreads)
  KNN_matrix = t(KNN_matrix) + KNN_matrix
  KNN_graph = graph_from_adjacency_matrix(KNN_matrix, weighted = TRUE, mode = "undirected")
  Clustering = cluster_louvain(KNN_graph, resolution = resolution)
  Clustering = membership(Clustering)
  Clustering = as.character(Clustering)
  return(Clustering)
}
