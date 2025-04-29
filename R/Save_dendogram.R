#' Save_dendogram
#'
#' This function generates and saves a dendrogram based on clustering results and Mean expression data.
#' It also creates a bar plot showing the frequency of annotations for each cluster and saves the 
#' combined output as a PDF.
#'
#' @param Clustering A vector or factor representing the clustering results.
#' @param Mean_expression A matrix or data frame containing mean expression values for each group.
#' @param Annotation_cells_renamed A vector of annotations for the cells.
#' @param Output_path A character string specifying the directory path where the output PDF plot will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#'
#' @return NULL This function does not return a value; it saves a PDF file.
#' @export
#' @import stats
#' @import dplyr
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' 
Save_dendogram = function (Clustering,Mean_expression,Annotation_cells_renamed,Output_path,Method,Tissue) {
  pdf(paste0(Output_path,'_',Method,'_',Tissue,"_Dendrogram.pdf"),width = 6.5,height = 6.5,useDingbats = FALSE)
  ratio <- table(Clustering, Annotation_cells_renamed)/rowSums(table(Clustering, Annotation_cells_renamed))
  filtered_table <- as.data.frame(ratio)
  filtered_table <- filtered_table %>%
    group_by(Clustering) %>%
    dplyr::filter(Freq == max(Freq))
  if (!is.data.frame(Expression_file)) {
    Mean_expression = t(Mean_expression)
  }
  correlation_matrix <- cor(Mean_expression, method = "spearman")
  distance_matrix <- as.dist(1 - correlation_matrix)
  hc <- hclust(distance_matrix)
  desired_order <- hc$labels[hc$order]
  filtered_table <- filtered_table[order(factor(filtered_table$Clustering, levels = desired_order)), ]
  par(las=1,mfrow=c(2,1))
  filtered_table <- filtered_table[order(factor(filtered_table$Clustering, levels = desired_order)), ]
  filtered_table$Clustering <- factor(filtered_table$Clustering, levels = unique(filtered_table$Clustering))
  m <- ggplot(filtered_table, aes(x = Clustering, y = Freq)) +
    geom_bar(stat = "identity")+
    ylim(0, 1) +
    scale_x_discrete(labels = filtered_table$Annotation_cells_renamed) +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  grid.arrange( plot(hc), m, ncol=1, heights=c(0.3, 0.5))
  dev.off()
  return(invisible())
}

