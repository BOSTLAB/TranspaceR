#' Convert numeric values to a color gradient from white to red.
#'
#' @param x A numeric vector to be converted to colors.
#' @param max_scale An optional maximum scale value for normalization.
#' @return A vector of colors corresponding to the input values.
#' @export
#' 

color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}
