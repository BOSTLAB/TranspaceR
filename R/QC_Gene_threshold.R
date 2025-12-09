#' Quality Control Gene Thresholding
#'
#' This function performs quality control on gene expression data by applying Otsu's thresholding method
#'  and uses quality probes to calculate a separability score.
#'
#' @param Expression_file A data frame containing gene expression data.
#' @param Meta_data A data frame containing metadata including radius information.
#' @param Method A character string specifying the method used for analysis.
#' @param Tissue A character string indicating the type of tissue being analyzed.
#' @param Output_path A character string specifying the directory where the output PDF will be saved.
#' @return A list containing the Otsu threshold and the separability score.
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @export

QC_Gene_threshold = function(Expression_file,Meta_data,Method,Tissue,Output_path) {
  
  radius_squared = Meta_data$radius^2
  Gene_size = Matrix::colSums(Expression_file)
  Lib_size = Matrix::rowSums(Expression_file)
  Neg_prob = (which(grepl(names(Gene_size),pattern = "Neg"))) #find positions of the negative probes 
  False_prob = (which(grepl(names(Gene_size),pattern = "False")))
  Blank = (which(grepl(names(Gene_size),pattern = "Blank")))
  Neg_prob = c(Neg_prob,False_prob,Blank)
  Neg_prob_count = log10(Gene_size[Neg_prob])
  Neg_prob_count = Neg_prob_count[is.finite(Neg_prob_count)]
  x = log10(Gene_size)
  x = x[is.finite(x)]
  threshold = Otsu_thresholding(x)
  neg_probes_before_threshold <- sum(Neg_prob_count < threshold)
  total_neg_probes <- length(Neg_prob_count)
  score_separability = neg_probes_before_threshold/total_neg_probes
  pdf(paste0(Output_path,Method,"_",Tissue,"_negprobes.pdf"),width = 6.5,height = 6.5,useDingbats = FALSE)
  par(fig=c(0,1,0,1),las=1,bty='l',mar=c(10,6,4,4),cex.lab = 0.8)
  y = hist(x,50,xaxs='i',yaxs='i',col = "lightblue",main=paste0(Method,"_",Tissue),
           cex.lab=1.5,xlab="",xaxs='i',yaxs='i',xlim=range(rbind(min(x),Neg_prob_count,max(x))))  
  abline(v = threshold, col = "red", lwd = 2)
  legend("topleft",legend = paste("Otsu's Threshold:", round(10^threshold)), col = "red", lwd = 2,xjust=1,yjust=1,cex = 0.7)
  legend("topright",legend = paste("Separability score:", round(score_separability,2)), col = "red", lwd = 2,xjust=1,yjust=1,cex = 0.7)
  mtext("Number of RNA molecules per gene (log10)", side = 1, line = 7, cex = 1.2, font = 2)
  par(fig=c(0,10,1.5,10)/10,bty='l',mar=c(5,6,4,4))
  par(fig=c(0,10,0,3.5)/10,bty='o')
  par(new=T)
  plot(NULL,xlim=range(rbind(min(x) ,Neg_prob_count,max(x))),xaxt='n',yaxt='n',ylim=c(0,1),xlab="",ylab="Neg probes",xaxs='i',yaxs='i')
  abline(v=Neg_prob_count,lty=1,col="red")
  dev.off()
  return(list(Otsu_threshold = 10^threshold,Separability_score=score_separability))
}

