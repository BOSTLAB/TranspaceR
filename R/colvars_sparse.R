#' Compute variance of columns of a sparse matrix
#'
#' @param mat A sparse matrix
#'
#' @return A vector of values corresponding to the variance of each column of the input.
#' @export
#' 

colvars_sparse <- function(mat) {
  means = colMeans(mat)
  squared_diff = mat
  for (j in 1:ncol(mat)) {
    start= mat@p[j]+1
    end= mat@p[j+1]
    if (start <= end) {
      idx=start:end
      squared_diff@x[idx]=(mat@x[idx]-means[j])^2 }}
  variances= numeric(ncol(mat))
  for (j in 1:ncol(mat)) {
    nnz_j= mat@p[j+1]- mat@p[j]
    n = nrow(mat)
    sum_sq_diff= sum(squared_diff@x[(mat@p[j]+1):mat@p[j+1]])
    mean_sq_diff=(sum_sq_diff+(n-nnz_j)*means[j]^2)/n
    variances[j]= mean_sq_diff }
  return(variances)
}
