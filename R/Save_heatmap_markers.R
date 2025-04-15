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
#' @import pheatmap
Save_heatmap_markers = function(Expression_file, object = Clustering,Output_path,Method,Tissue) {
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
  Log2FC_list = sapply(genes, function(gene) calculate_log2fc(gene,object), simplify = FALSE)
  Log2FC = as.data.frame(Log2FC_list)
  result <- Log2FC %>%
    rownames_to_column(var = "Cluster") %>%
    pivot_longer(-Cluster, names_to = "Gene", values_to = "logFC") %>%
    group_by(Cluster) %>%
    top_n(5, logFC) %>%
    ungroup() %>%
    pivot_wider(names_from = Gene, values_from = logFC)
  col = unique(colnames(result))[-1]
  Log2FC= Log2FC[,col]
  Log2FC = as.data.frame(Log2FC)
  pdf(paste(Output_path,Method,'_',Tissue,"_markers_heatmap.pdf",sep=''), useDingbats = FALSE, width = 12, height = 12)
  p <- pheatmap(t(scale(mean_expression[rownames(Log2FC),colnames(Log2FC)])), clustering_method = 'ward.D')
  print(p)
  dev.off()
  return(invisible())
}
