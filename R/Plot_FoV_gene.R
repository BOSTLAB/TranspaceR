#' Plot Field of View for a Gene
#'
#' This function plots the cell centroids colored by the expression of a specified gene.
#'
#' @param gene A string representing the gene for which the expression is to be plotted.
#' @return A plot of cell centroids colored by gene expression.
#' @export
#' 
Plot_FoV_gene = function(gene) {
  plot(Meta_data$cell_centroid_x,Meta_data$cell_centroid_y,
       pch=21,bg=color_convertion(Expression_file[,gene]),lwd=0.1,cex=0.1)
}
