#' Pad Definitor
#'
#' This function calculates an optimal padding size for rasterization based on the cell centroid coordinates in the metadata.
#'
#' @param Meta_data A data frame containing metadata with cell centroid coordinates.
#' @return An integer representing the calculated padding size.
#' @export
#' 
pad_definitor <- function (Meta_data) {
  xrange = max(Meta_data$cell_centroid_x)-min(Meta_data$cell_centroid_x) 
  yrange = max(Meta_data$cell_centroid_y)-min(Meta_data$cell_centroid_y) 
  range = mean(c(xrange,yrange))
  n_pad = round(range/200)
  if (n_pad == 0) { n_pad = 1}
  return(n_pad)
}
