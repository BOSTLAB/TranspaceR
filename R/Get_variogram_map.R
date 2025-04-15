#' Get Variogram Map
#'
#' This function computes the variogram map of a given image using the Fast Fourier Transform (FFT).
#'
#' @param im A matrix representing the image for which the variogram map is to be computed.
#' @return A matrix representing the variogram map of the input image.
#' @export
#' @import imager
#' @import gsignal
#' @import OpenImageR
Get_variogram_map = function(im) {
  im_index = im
  im_index[!is.na(im_index)] = 1
  im_index[is.na(im_index)] = 0
  
  im[is.na(im)]=0
  
  im_ft = fft(im)
  im_index_ft = fft(im_index)
  im_2_ft = fft(im^2)
  N = nrow(im)*ncol(im)
  
  Variogram_map = Re(fft(im_2_ft*Conj(im_index_ft)+Conj(im_2_ft)*im_index_ft-2*im_ft*Conj(im_ft),inverse = TRUE))/(2*Re(fft(im_index_ft*Conj(im_index_ft),inverse = TRUE)))
  
  Variogram_map_reshaped = fftshift(Variogram_map)
  return(Variogram_map_reshaped)
}

