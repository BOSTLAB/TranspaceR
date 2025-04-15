#' Plot Field of View with Clusters
#'
#' This function plots the cell centroids colored by their classification.
#' @param object Default is KNN clustering.
#' @return A plot of cell centroids colored based on the classification.
#' @export
#'
Plot_FoV_cluster = function(object = Clustering) {
  plot(Meta_data$cell_centroid_x,Meta_data$cell_centroid_y,
       pch=21,bg=string.to.colors(object),lwd=0.1,cex=0.2)
}

