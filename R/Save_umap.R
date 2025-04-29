#' Save UMAP Plot
#'
#' This function saves a UMAP plot of colored points, based on clustering results, to a PNG file.
#'
#' @param Umap_result A matrix containing UMAP results.
#' @param object A vector of clustering results corresponding to the UMAP points.
#' @param Output_path Character string specifying the output directory where to save the plot.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#'
#' @return NULL This function does not return a value; it saves a PNG file.
#' @importFrom uwot umap
#' @export

Save_umap <- function(Umap_result,object = Clustering,Output_path,Method,Tissue) {
  optimal_palette <- rainbow(length(unique(object)))
  names(optimal_palette) <- unique(object)
  png(paste0(Output_path, Method, "_", Tissue, "_umap.png"), width = 1000, height = 1000)
  plot(Umap_result, pch = 21, bg = optimal_palette[object], cex = 1, lwd = 0.2)
  legend("topleft", legend = unique(object), fill = optimal_palette[unique(object)], lty = 1, cex = 1, box.lty = 0.4, box.lwd = 0.4)
  dev.off()
  return(invisible())
}
