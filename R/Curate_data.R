#' Curate expression data and metadata based on quality control thresholds.
#'
#' @param Expression_file A matrix or dataframe of  gene expression data.
#' @param Meta_data A data frame containing metadata for each cell.
#' @param QC2_results A data frame containing quality control results.
#' @param min_lib_size Minimum total number of transcripts per cell threshold.
#' @param max_lib_size Maximum total number of transcripts per cell threshold.
#' @param min_cell_radius Minimum cell radius threshold.
#' @param max_cell_radius Maximum cell radius threshold.
#' @return A list containing curated expression data and metadata.
#' @export
#' 
Curate_data = function (Expression_file,Meta_data,QC2_results,min_lib_size,max_lib_size,min_cell_radius,max_cell_radius) {
  Gene_size = colSums(Expression_file)
  Lib_size = rowSums(Expression_file)
  radius = Meta_data$radius
  threshold = QC2_results$Otsu_threshold
  Selected_cells  = Lib_size>min_lib_size & Lib_size<max_lib_size & radius >min_cell_radius & radius < max_cell_radius
  Selected_genes = Gene_size > threshold
  Expression_file = Expression_file[Selected_cells,Selected_genes]
  Meta_data = Meta_data[Selected_cells,]
  return(list(Expression_file = Expression_file, Meta_data = Meta_data))
}
