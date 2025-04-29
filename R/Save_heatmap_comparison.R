#' Save_heatmap_comparison
#'
#' This function generates and saves a heatmap comparing clustering results with cell annotations.
#' It calculates the frequency of annotations for each cluster and visualizes the results in a heatmap.
#'
#' @param Clustering A vector representing the clustering results.
#' @param Annotation_cells_renamed A vector of annotations for the cells.
#' @param Output_path A character string specifying the directory where the output PDF will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#'
#' @return  NULL This function does not return a value; it saves a PDF file.
#' @export
#' @importFrom pheatmap pheatmap 

Save_heatmap_comparison = function (Clustering,Annotation_cells_renamed,Output_path,Method,Tissue) {
  ratio <- table(Clustering, Annotation_cells_renamed)/rowSums(table(Clustering, Annotation_cells_renamed))
  filtered_table <- as.data.frame(ratio)
  filtered_table <- filtered_table %>%
    group_by(Clustering) %>%
    dplyr::filter(Freq == max(Freq))
  table_structure <- xtabs(Freq ~ Clustering + Annotation_cells_renamed, data = filtered_table)
  # KNN vs Scimilarity
  pdf(paste(Output_path, Method,'_',Tissue, "_Heatmap_KNN_SCI.pdf", sep=""), useDingbats = FALSE, width = 8, height = 8)
  p <- pheatmap(table_structure, clustering_method = 'ward.D', color = colorRampPalette(c("beige", "orange", "red"))(100),  breaks = seq(0, 1, length.out = 100), fontsize = 10)
  print(p)
  dev.off()
  return(invisible())
}
