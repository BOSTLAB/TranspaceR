#' Select and visualize shared genes
#'
#' This function takes a list of gene sets and generates an UpSet plot to visualize the shared genes among them.
#'
#' @param Selected_objects A list containing vectors of gene names. Each vector represents a set of genes (e.g., variance genes, zero genes).
#' @param Selected_names A character vector containing names for the gene sets corresponding to the vectors in `Selected_objects`.
#'
#' @return A vector of unique shared genes across the provided gene sets.
#'
#' @import UpSetR
#' @import grDevices
#' @import utils
#' @export

Select_genes = function (Selected_objects =list(Variance_genes,Zero_genes),Selected_names = c('Variance_genes','Zero_genes')) {
  Shared_genes = unique(unlist(Selected_objects))
  #Upset plot
  score_df = as.data.frame(Shared_genes)
  for (i in 1:length(Selected_names)) {
    score_df[[Selected_names[i]]] <- as.integer(score_df$Shared_genes %in% Selected_objects[[i]])
  }
  pdf(paste0(Output_path, Method, "_", Tissue, "_Upset-plot.pdf"), width = 6.5, height = 6.5, useDingbats = FALSE)
  p <- upset(score_df, sets = Selected_names,
             order.by = "freq", empty.intersections = "on", text.scale = 1.5, matrix.color = "#56B4E9",
             main.bar.color = "grey", mainbar.y.label = 'Number of shared genes')
  print(p)
  dev.off()
  return(Shared_genes)
}