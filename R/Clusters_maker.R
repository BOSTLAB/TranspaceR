#' Create clusters and compute log2 fold changes for shared genes in expression data.
#'
#' @param Expression_file A data frame or matrix of gene expression data.
#' @param Shared_genes A vector of gene names that are shared across samples.
#' @param K An integer specifying the number of neighbors for clustering (default is 30).
#' @param metric_used A string indicating the distance metric to use (default is "L2").
#' @param nThreads An integer specifying the number of threads for computation (default is 12).
#' @param resolution A numeric value for the resolution parameter in clustering (default is 1).
#' @return A list containing PCA data, corrected data, cluster assignments, mean expression, and log2 fold change table.
#' @export
#' @import dplyr
#' @import irlba
#' @import tibble
Clusters_maker = function(Expression_file, Shared_genes, K = 30, metric_used = "L2", nThreads = 12, resolution = 1) {
  if (is.data.frame(Expression_file)){
    Data_correction = Expression_file[, Shared_genes]
    Data_correction = as.data.frame(Data_correction)
    PCA_data = irlba(as.matrix(Data_correction), nv = 50, verbose = TRUE)
    Clustering = Cell_clustering_function(PCA_data$u, K, metric_used, nThreads, resolution)
    Mean_expression = aggregate(Data_correction, by = list(Clustering), FUN = mean)
    rownames(Mean_expression) = Mean_expression$Group.1
    Mean_expression = t(Mean_expression[,-1])
    genes = rownames(Mean_expression)
  } else {
    Data_correction = Expression_file[, Shared_genes]
    PCA_data = irlba(as.matrix(Data_correction), nv = 50, verbose = TRUE)
    Clustering = Cell_clustering_function(PCA_data$u, K, metric_used, nThreads, resolution)
    Mean_expression = aggregate_sparse(Data_correction,Clustering)
    genes = rownames(Mean_expression)
    Mean_expression = t(Mean_expression[,-1])
  }
  Log2FC_list = sapply(genes, function(gene) calculate_log2fc(gene, Clustering), simplify = FALSE)
  Log2FC = as.data.frame(Log2FC_list)
  result = data.frame(matrix(ncol = 3, nrow = ncol(Log2FC)*nrow(Log2FC)))
  colnames(result) = c('Cluster','Gene','LogFC')
  result$Cluster = rownames(Log2FC)
  result$LogFC = unlist(Log2FC)
  result$Gene = rep(colnames(Log2FC), each = nrow(Log2FC))
  result <- result %>%
    group_by(Cluster) %>%
    top_n(5, LogFC) %>%
    ungroup()
  result_wide = data.frame(matrix(ncol = length(unique(result$Gene))+1, nrow = length(unique(result$Cluster))))
  colnames(result_wide) = c('Cluster',unique(result$Gene))
  result_wide$Cluster = unique(result$Cluster)
  logfc_matrix <- tapply(result$LogFC, list(result$Cluster, result$Gene), FUN = function(x) x)
  result = as.data.frame(logfc_matrix)
  result = data.frame(Cluster = rownames(logfc_matrix),result)
  col = unique(colnames(result))[-1]
  Log2FC = Log2FC[, col]
  Log2FC = as.data.frame(Log2FC)
  return(list(PCA_data = PCA_data, Data_correction = Data_correction, Clustering = Clustering, Mean_expression = Mean_expression, Log2FC_table = Log2FC))
}
