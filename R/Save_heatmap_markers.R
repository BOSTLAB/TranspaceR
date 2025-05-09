#' Save Heatmap Markers
#'
#' This function generates a heatmap of the top marker genes for each cluster based on the provided expression data.
#' It calculates the log2 fold change for each gene and selects the top 5 genes with the highest log2 fold change
#' for each cluster. The resulting heatmap is saved as a PDF file.
#'
#' @param Expression_file A data frame or a sparse matrix containing gene expression data with genes as columns and cells as rows.
#' @param object A clustering object that defines the clusters for the samples. This is typically a vector 
#' or factor indicating the cluster assignment for each sample.
#' @param Output_path A character string specifying the directory where the output PDF file will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#' @return  NULL This function does not return a value; it saves a PDF file.
#' @export
#' @importFrom pheatmap pheatmap 
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
Save_heatmap_markers = function(Expression_file, object = Clustering, Output_path,Method,Tissue,name_object = 'Clustering') {
  if (is.data.frame(Expression_file)) {
  mean_expression = aggregate(Expression_file, FUN = mean, by = list(object))
  rownames(mean_expression) <- mean_expression[[1]]
  mean_expression <- mean_expression[-1]
  genes = colnames(mean_expression)
  }
  else{
    mean_expression = aggregate_sparse(Data_correction,object)
    genes = rownames(mean_expression)
    mean_expression = t(mean_expression)
  }
  Log2FC_list = sapply(genes, function(gene) calculate_log2fc(gene, object), simplify = FALSE)
  Log2FC = as.data.frame(Log2FC_list)
  result = data.frame(matrix(ncol = 3, nrow = ncol(Log2FC)*nrow(Log2FC)))
  colnames(result) = c('Cluster','Gene','LogFC')
  result$Cluster = rownames(Log2FC)
  result$LogFC = unlist(Log2FC)
  result$Gene = rep(colnames(Log2FC), each = nrow(Log2FC))
  result <- result %>%
    group_by(Cluster) %>%
    top_n(5, LogFC) %>%
    ungroup()
  result_wide = data.frame(matrix(ncol = length(unique(result$Gene))+1, nrow = length(unique(result$Cluster))))
  colnames(result_wide) = c('Cluster',unique(result$Gene))
  result_wide$Cluster = unique(result$Cluster)
  logfc_matrix <- tapply(result$LogFC, list(result$Cluster, result$Gene), FUN = function(x) x)
  result = as.data.frame(logfc_matrix)
  result = data.frame(Cluster = rownames(logfc_matrix),result)
  col = unique(colnames(result))[-1]
  Log2FC = Log2FC[, col]
  Log2FC = as.data.frame(Log2FC)
  pdf(paste0(Output_path,Method,'_',Tissue,'_', name_object, "_markers_heatmap.pdf"), useDingbats = FALSE, width = 12, height = 12)
  p <- pheatmap(t(scale(mean_expression[rownames(Log2FC),colnames(Log2FC)])), clustering_method = 'ward.D')
  print(p)
  dev.off()
  return(invisible())
}
