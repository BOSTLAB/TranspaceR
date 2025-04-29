#' Save Boxplot of Gene Expression
#'
#' This function generates and saves a boxplot of gene expression across different cell types.
#'
#' @param Data_correction A data frame containing corrected gene expression data.
#' @param object A vector of cell type assignments.
#' @param gene A string representing the gene for which the boxplot is to be created.
#' @param Output_path A character string specifying the directory path where the output PDF plot will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#' @return  NULL This function does not return a value; it saves  a PDF file.
#' @export
#' 
 

Save_boxplot = function(Data_correction,object = Annotation_cells_renamed,gene ='MS4A1',Output_path,Method,Tissue) {
  pdf(paste0(Output_path,Method,'_',Tissue,"_expression_boxplot.pdf"),useDingbats = FALSE,width = 6.5,height = 6.5)
  par(las=1,bty="l",mar=c(5,15,4,3))
    boxplot(Data_correction[,gene]~object,
            outline=F,horizontal = TRUE,ylab="",xlab="Gene expression",cex.lab=1.3,main=gene)
  dev.off()
  return(invisible())
}

