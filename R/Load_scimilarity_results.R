#' Load Scimilarity Results
#'
#' This function loads cell type predictions from a CSV file and normalizes the counts.
#' Cell types below a specified threshold are grouped into an "Other" category.
#'
#' @param path A string representing the file path to the directory containing the CSV file.
#' @param celltype_threshold A numeric value representing the threshold for cell type normalization.
#' @return A vector of renamed cell types.
#' @export
#' 
Load_scimilarity_results= function(path,celltype_threshold=0.01) {
  Annotation_file = read.delim(paste(path,"cell_predictions.csv",sep = ""), sep=",")
  Annotation_cells = Annotation_file$predictions_unconstrained
  count_cell = table(Annotation_cells)
  count_cell_normalised = count_cell/sum(count_cell)
  celltype_to_remove = names(which(count_cell_normalised<celltype_threshold))
  Annotation_cells_renamed = Annotation_cells
  Annotation_cells_renamed[Annotation_cells_renamed%in%celltype_to_remove] = "Other"
  return(Annotation_cells_renamed)
}
