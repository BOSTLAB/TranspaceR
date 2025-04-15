#' Calculate log2 fold change for a given gene across clusters.
#'
#' @param gene A string representing the gene of interest.
#' @param Clustering A vector indicating the cluster membership of each sample.
#' @return A named vector of log2 fold changes for the specified gene.
#' @export
#' 
calculate_log2fc <- function(gene,Clustering) {
  x = Expression_file[, gene]
  mean_expression = aggregate(x, FUN = mean, by = list(Clustering))$x
  Proportion_cluster = table(Clustering) / length(Clustering)
  
  Matrix_proportion = matrix(Proportion_cluster, ncol = length(Proportion_cluster), nrow = length(Proportion_cluster), byrow = TRUE)
  diag(Matrix_proportion) = 0
  Matrix_proportion = Matrix_proportion / rowSums(Matrix_proportion)
  
  Mean_expression_other_cluster = as.numeric(mean_expression %*% t(Matrix_proportion))
  Log2FC = log2(mean_expression / Mean_expression_other_cluster)
  names(Log2FC) = names(Proportion_cluster)
  
  return(Log2FC)
}

