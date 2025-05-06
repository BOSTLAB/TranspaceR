#' Compute variance of columns of a sparse matrix
#'
#' @param mat A sparse matrix
#'
#' @return A vector of values corresponding to the variance of each column of the input.
#' @export
#' 

colvars_sparse <- function(mat) {
  variances = apply(mat, 2, function(col_data) {
    col_mean = mean(col_data, na.rm = TRUE)
    mean((col_data - col_mean)^2, na.rm = TRUE)
  })
  return(variances)
}
