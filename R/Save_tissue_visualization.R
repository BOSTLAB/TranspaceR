#' Save Tissue Visualization
#'
#' This function generates a visualization of tissue based on cell centroids and clusters.
#' It saves the plot as a PNG file in the specified output path.
#'
#' @param Meta_data A data frame containing the metadata for the cells with
#' the x and y coordinates (columns named `cell_centroid_x` and `cell_centroid_y`).
#' @param object A vector representing classification of cells. Default is KNN Clustering result.
#' @param Output_path A character string specifying the directory where the output PNG file will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#' @param scaling_factor A numeric value that determines the size of the points in the plot. Default is 1.5.
#'
#' @return NULL This function does not return a value; it saves a PNG file.
#' @export
#' @import spatstat 
Save_tissue_visualization = function(Meta_data, object = Clustering, Output_path,Method,Tissue,scaling_factor=1.5) {
  ppp_temp = ppp(x = Meta_data$cell_centroid_x,y = Meta_data$cell_centroid_y,
                 window = owin(xrange = range(Meta_data$cell_centroid_x),yrange = range(Meta_data$cell_centroid_y)))
  N_cluster = length(unique(object))
  optimal_palette = rainbow(N_cluster)
  png(paste0(Output_path,Method,"_",Tissue,"_atlas.png"),width = 4000,height = 4000)
  plot(ppp_temp[1,],main="",bty="n",cex = 0)
  for (k in 1:N_cluster) {
    temp_cell_type = unique(object)[k]
    points(ppp_temp[object==temp_cell_type,],cex=scaling_factor,col=optimal_palette[k])
  }
  dev.off()
}
