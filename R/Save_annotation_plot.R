#' Save Annotation Plot
#'
#' This function generates a bar plot showing the frequency of different annotation categories
#' in the provided data and saves it as a PDF file.
#'
#' @param Annotation_cells_renamed A vector of annotation categories for cells. This should contain the annotations to visualize.
#' @param Output_path A character string specifying the directory path where the output PDF file will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#' @return  NULL This function does not return a value; it saves a PDF file.
#' @import ggplot2
#' @import ggthemes
#' @export
#' 
#' 
Save_annotation_plot = function (Annotation_cells_renamed,Output_path,Method,Tissue){
  category_freq= (table(Annotation_cells_renamed)/length(Annotation_cells_renamed))*100
  category_freq = as.data.frame(category_freq)
  colnames(category_freq) = c('Annotation_cells','Freq')
  category_freq$Freq = round(category_freq$Freq,2)
  category_freq$Annotation_cells <- factor(category_freq$Annotation_cells, levels = category_freq$Annotation_cells[order(category_freq$Freq, decreasing = FALSE)])
  
  pdf(paste0(Output_path,Method,"_",Tissue,"_annotation_proportion.pdf"),width = 6.5,height = 6.5,useDingbats = FALSE)
  par(fig=c(0,1,0,1),las=1,bty='l',cex.lab = 0.7, cex.axis = 0.7)
  p <- ggplot(data = category_freq, aes(x = Annotation_cells, y = Freq)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = paste0(Method,'_',Tissue)) +
    theme_clean() +
    theme(axis.text = element_text(size = 10, angle = 0, hjust = 1),
          axis.title = element_text(size = 14)) +
    coord_flip()+
    geom_text(aes(label=Freq), vjust=0.5, size=3)
  print(p)
  dev.off()
  return(invisible())
}
