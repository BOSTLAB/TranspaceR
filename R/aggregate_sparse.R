#' aggregate_sparse
#'
#' Aggregate a sparse matrix to compute the mean for each category of an object
#'
#' @param Expression_file Gene Expression file as sparse matrix,  
#' @param object A vector containing clustering or annotation
#' @return Mean expression of genes in each group matrix
#' @export
#' @import Matrix
aggregate_sparse <- function(Expression_file, object = Clustering) {
  groups=unique(object)
  genes =colnames(Expression_file)
  mean_expression= matrix(NA, nrow = length(genes), ncol = length(groups))
  rownames(mean_expression)= genes
  colnames(mean_expression) = as.character(groups)
  objectr = as.factor(object)
  mean_expression = sapply(groups, function(g) {
    group_indices = which(object == g)
    if (length(group_indices) > 0) {
      colMeans(Expression_file[group_indices, , drop = FALSE], na.rm = TRUE)
    } else {
      NA 
    }
  })
  return(mean_expression)
}
