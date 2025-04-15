#' Plot Field of View
#'
#' This function plots the cell centroids for a specified field of view (FoV).
#'
#' @param FoV An integer representing the field of view to be plotted.
#' @return A plot of cell centroids for the specified FoV.
#' @export
#' 
Plot_FoV = function(FoV = 1) {
  plot(Meta_data$cell_centroid_x[Meta_data$dataset.suffix==FoV],Meta_data$cell_centroid_y[Meta_data$dataset.suffix==FoV],
       pch=21,lwd=0.1,cex=0.1,bg=string.to.colors(Clustering[Meta_data$dataset.suffix==FoV]))
}

