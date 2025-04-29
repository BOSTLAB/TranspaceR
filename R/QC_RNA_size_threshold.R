#' QC_RNA_size_threshold
#'
#' This function performs quality control on RNA expression data by analyzing the relationship between total RNA and cell size (radius).
#' It generates a 2D plot of RNA expression versus cell size and saves it as a PDF.
#'
#' @param Expression_file A matrix or data frame containing gene expression data.
#' @param Meta_data A data frame containing metadata including a column named 'radius'.
#' @param Method A character string specifying the method used for analysis.
#' @param Tissue A character string indicating the type of tissue being analyzed.
#' @param Output_path A character string specifying the directory where the output PDF will be saved.
#' @return  NULL This function does not return a value; it saves a PDF file.
#' @export
#' @import ggplot2
#' @import ggthemes
#' 
QC_RNA_size_threshold = function(Expression_file,Meta_data,Method,Tissue,Output_path) {
  
  radius_squared = Meta_data$radius^2
  Lib_size = rowSums(Expression_file)
  # Saving QC plots 
  # 2D plot of RNA~Size
  m = stats::lm((Lib_size)~ 0 + radius_squared,subset = Lib_size>0)
  x_range = Meta_data$radius
  y_fitted = coef(m)*radius_squared
  m_0 = summary(m)
  m_0$r.squared
  # Create dataframe for visualization
  Lib_df = data.frame(radius = Meta_data$radius,Lib_size = Lib_size)
  # Visualization
  pdf(paste0(Output_path,Method,"_",Tissue,"_2D.pdf"),width = 7,height = 7,useDingbats = FALSE)
  p <- ggplot(data = Lib_df, aes(x = radius, y = Lib_size)) +
    geom_bin2d(bins = 150) +
    scale_y_continuous(trans = "log",breaks=seq(0,max(Lib_size),round(max(Lib_size)/10,-2)+50)) + 
    scale_fill_continuous(type = "viridis") +
    theme_bw()+
    xlab("Cell radius (um)") +
    ylab("Total RNA") +
    labs(title = paste0(Method,"_",Tissue)) +
    theme_clean() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10))+
    geom_line(data = data.frame(x = x_range, y = y_fitted), aes(x = x, y = y), color = "blue", linetype = "dashed", size = 1.5)+
    annotate("text",x=Inf,y=Inf, label = paste("R-squared = ", round(m_0$r.squared,2), sep =""), size = 4, color = "blue",hjust=1,vjust=1)
  print(p)
  dev.off()
  return(invisible())
}
