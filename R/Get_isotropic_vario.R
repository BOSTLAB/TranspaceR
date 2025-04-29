#' Get the isotropic variogram from a given variogram matrix.
#'
#' @param Variogram A matrix representing the variogram.
#' @return A data frame with distances and corresponding variogram values.
#' @export

Get_isotropic_vario = function(Variogram) {
  n_x = ncol(Variogram)
  n_y = nrow(Variogram)
  
  center_x = round(n_x/2)
  center_y = round(n_y/2)
  
  m_x = t(apply(matrix(0,ncol = n_x,nrow = n_y),MARGIN = 1,FUN = function(x) {-center_x:(center_x-1)}))
  m_y = apply(matrix(0,ncol = n_x,nrow = n_y),MARGIN = 2,FUN = function(x) {-center_y:(center_y-1)})
  
  adjust_dimensions <- function(m, dim_var) {
    dim= dim(m)
    
    if (dim[1] != dim_var[1]) {
      m_adjusted = m[rep(1:dim[1], length.out = dim_var[1]), ]
    } else {
      m_adjusted = m  
    }
    if (dim[2] != dim_var[2]) {
      m_adjusted= m_adjusted[, rep(1:dim[2], length.out = dim_var[2])]
    }
    return(m_adjusted)
  }
  
  dim_var <- dim(Variogram)
  
  m_x= adjust_dimensions(m_x, dim_var)
  m_y= adjust_dimensions(m_y, dim_var)
  
  m_dist = sqrt(m_x^2+m_y^2)
  u = aggregate(as.numeric(Variogram),by = list(as.numeric(m_dist)),FUN = mean)
  colnames(u) = c("r","Variogram")
  return(u)
  
}
