#' Otsu Thresholding
#'
#' This function implements in-house Otsu's method to determine an optimal threshold for a given dataset.
#'
#' @param x A numeric vector of data for which the threshold is to be determined.
#' @param number_bins An integer specifying the number of bins to use for histogram calculation.
#' @return A numeric value representing the selected threshold.
#' @export
#' 
Otsu_thresholding = function(x,number_bins=100) {
  List_bin = quantile(x,base::seq(from=0,to=1,length.out=number_bins))
  Intravariance_vector = c()
  for (k in 1:number_bins) {
    threshold_temp = List_bin[k]
    s = length(x[x<threshold_temp])*var(x[x<threshold_temp]) + length(x[x>threshold_temp])*var(x[x>threshold_temp])
    Intravariance_vector= c(Intravariance_vector,s)
  }
  Selected_values = List_bin[which.min(Intravariance_vector)]
  return(Selected_values)
}
