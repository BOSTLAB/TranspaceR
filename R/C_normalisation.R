#' Normalize cell data based on cell size.
#'
#' @param y A matrix of expression data.
#' @param scaling_factor A vector representing the area of each cell.
#' @return A matrix of normalized cell data.
#' @export
#' 
C_normalisation = function(y, scaling_factor){
  Cell_size = scaling_factor  # Area of cell
  Tau_parameter = colMeans(y/Cell_size)
  Scaling_parameter = matrix(Cell_size,ncol = ncol(y),nrow = nrow(y),byrow = FALSE)
  Scaling_parameter = t(apply(Scaling_parameter,MARGIN = 1,FUN = function(x) {x + 1/Tau_parameter}))
  Scaled_data = y / Scaling_parameter
  return(Scaled_data)
}
