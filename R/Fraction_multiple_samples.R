#' This function takes the expression file and metadata and separates them into multiple ss based on the
#' metadata provided. It returns a list containing the expression data and metadata for each unique s.
#' @param Expression_file 
#' @param Meta_data It must include a column named 's' that identifies the ss corresponding to the expression data.
#'
#' @return A list containing two elements:
#' \item{Expression_files}{A list of data frames, each corresponding to a unique s, filtered from the Expression_file.}
#' \item{Meta_data_files}{A list of data frames, each corresponding to the metadata for each unique s.}
#' @export

Fraction_multiple_samples <- function(Expression_file, Meta_data) {
  if (!"sample" %in% colnames(Meta_data)) {
    stop("No 's' column provided in Metadata")
  }
  expression_files = list()
  meta_data_files = list()
  samples=unique(Meta_data$sample)
  for (s in samples) {
    idx <- which(Meta_data$sample == s)
    expression_files[[s]] = Expression_file[idx, , drop = FALSE]
    meta_data_files[[s]] = Meta_data[idx, , drop = FALSE]
  }
  return(list(Expression_files = expression_files, Meta_data_files = meta_data_files))
}

