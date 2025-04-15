#' Save Geary Variance Plot
#'
#' This function generates a scatter plot comparing Variance scores and Geary scores for Selected genes.
#' The plot highlights Shared genes, Variance-only genes, and Geary-only genes with different colors.
#'
#' @param Variance_computation A list containing the results of variance computation.
#' @param Geary_computation A list containing the results of Geary computation. 
#' @param Output_path A character string specifying the directory path where the output PDF plot will be saved.
#' @param Method A character string indicating the spatial method used.
#' @param Tissue A character string indicating the type of tissue.
#' @return NULL This function does not return a value; it saves a PDF file.
#' @export
#' 

Save_geary_variance_plot = function(Variance_computation,Geary_computation,Output_path,Method,Tissue) {
  Variance_genes= Variance_computation$Selected_genes
  Variance_score = Variance_computation$Excess_variance_ratio
  Geary_score = Geary_computation$Values
  Geary_genes = Geary_computation$Selected_genes
  Shared_genes=intersect(Variance_genes,Geary_genes)
  Geary_genes_only = Geary_genes[!Geary_genes %in% Shared_genes]
  Variance_genes_only = Variance_genes[!Variance_genes %in% Shared_genes]
  pdf(paste0(Output_path,Method,"_",Tissue,"_spatial-variance.pdf"),width = 6.5,height = 6.5,useDingbats = FALSE)
  par(fig=c(0,1,0,1),las=1,bty='l',cex.lab = 0.7, cex.axis = 0.7)
  plot(Variance_score, Geary_score,xlab="Variance Score",ylab="Spatial Score",pch=21,col='lightgray')
  points(Variance_score[Shared_genes], Geary_score[Shared_genes],pch=21, col='black' ,bg='red', lwd = 1)
  points(Variance_score[Variance_genes_only], Geary_score[Variance_genes_only],pch=21,col='black', bg='orange', lwd = 1)
  points(Variance_score[Geary_genes_only], Geary_score[Geary_genes_only],pch=21, col='black',bg='green', lwd = 1)
  
  text(Variance_score[Variance_genes][1:20], Geary_score[Variance_genes][1:20],labels=Variance_genes[1:20], pos=4, cex=0.4, col='black')
  text(Variance_score[Geary_genes][1:20], Geary_score[Geary_genes][1:20],labels=Geary_genes[1:20], pos=4, cex=0.4, col='black')
  legend("topright",legend=c('Shared Genes', 'Variance Genes', 'Geary Genes'),col=c('red', 'orange', 'green'),pch=21,pt.bg=c('red', 'orange', 'green'),bty='n')
  m_0 <- lm(Geary_score ~ Variance_score)
  abline(m_0$coefficients[1], m_0$coefficients[2], col = "red", lty = 2, lwd = 2)
  model_sum = summary(m_0)
  legend("bottomright", legend = round(model_sum$r.squared, 3), col = "red", lwd = 2, xjust = 1, yjust = 1, cex = 0.6)
  title(main = paste0(Method,"_",Tissue), cex.sub = 0.6)
  dev.off()
  return(invisible())
}

