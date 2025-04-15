#' Excess Zero Score computation using Negative Binomial Distribution
#'
#' This function analyzes gene expression data to identify genes exhibiting an excess of zero counts.
#' It employs a negative binomial model to calculate p-values and selects significant genes based on
#' user-defined thresholds for p-values and delta scores and generates a plot of the results.
#'
#' @param Expression_file A data frame or matrix containing gene expression data, where rows represent genes
#' and columns represent samples.
#' @param Method A character string specifying the method used for analysis.
#' @param Tissue A character string indicating the type of tissue being analyzed.
#' @param Output_path A character string specifying the directory where the output PDF will be saved.
#' @param P_value_threshold A numeric value representing the threshold for statistical significance (default is 0.01).
#' @param Delta_threshold A numeric value representing the threshold for the excess zero score (default is 0.01).
#'
#' @return A list containing:
#'    - Selected_genes: A data frame of genes that show a significant excess of zeros based on the specified thresholds.
#'    - Excess_zero_score: A numeric vector of excess zero scores for each gene in the expression dataset.
#'
#' @export
#' @import MASS
#' @import N2R
#' 
Excess_zero_score_NB = function(Expression_file,Output_path,Method,Tissue,P_value_threshold = 0.01,Delta_threshold = 0.01) {
  
  Total_gene_expression = colSums(Expression_file)
  Proportion_zero = colSums(Expression_file==0)/nrow(Expression_file)
  Mean_expression = colMeans(Expression_file)
  Mean_expression[is.na(Mean_expression)]
  Var_expression = apply(Expression_file,FUN = var,MARGIN = 2)
  #Basic model using a NB distribution as background hypothesis (with a constant theta/over-dispersion parameter)
  #The probability of observing a zero is (theta/(mu+theta))^theta, for simpler fitting : log transform
  
  theta_init = median(1/((Var_expression-Mean_expression)/(Mean_expression^2)))
  m_null = nls(log(Proportion_zero) ~ theta*log(theta)-theta*log(theta+Mean_expression),algorithm = "port",
               lower =0,start = list(theta=theta_init),subset = Proportion_zero>0 & Proportion_zero<1) 
  theta = coef(m_null)
  
  Expected_proportion_zeros = (theta/(theta+Mean_expression))^theta
  Delta_excess_zero = Proportion_zero-Expected_proportion_zeros
  Mu_values=quantile(Mean_expression,probs = seq(0,1,length.out = 30))
  
  N_cells = nrow(Expression_file)
  N_cell_threshold = 10000
  Too_many_cells = FALSE
  
  N_simulations = 500
  Estimated_s = coef(m_null)
  
  if (N_cells>N_cell_threshold) {
    Too_many_cells = TRUE
  }
  
  if (!Too_many_cells) {
    cat("Computing the regular confidence interval....")
    Table_variance = c()
    
    for (i in 1:length(Mu_values)) {
      x = matrix(1,nrow = 1,ncol = N_simulations)
      x = apply(x,MARGIN = 2,FUN = function(x) {rnbinom(mu = Mu_values[i],n =N_cells,size = 1/Estimated_s )})
      Variance_estimated = apply(x,MARGIN = 2,FUN = var)
      Variance_estimator  = var(Variance_estimated) #Variance of the variance estimator : don't get lost ! 
      Table_variance = rbind(Table_variance,Variance_estimator)
    }
    rownames(Table_variance) = c()
    cat(" done ! \n")
    
  }
  if (Too_many_cells) {
    cat("Too many cells... confidence interval computed by sub-sampling")
    
    Table_variance <- foreach(i = 1:length(Mu_values), .combine = 'rbind') %dopar% {
      
      Variance_estimated = c()
      N_cells_values = c(1000,1500,2000,2500,3000,5000,10000)
      
      Table_variance_temp <- c()
      for (j in 1:length(N_cells_values)) {
        x = matrix(1,nrow = 1,ncol = N_simulations)
        x = apply(x,MARGIN = 2,FUN = function(x) {rnbinom(mu = Mu_values[i],n =N_cells_values[j],size = 1/Estimated_s)})
        Proportion_zero_temps = apply(x,MARGIN = 2,FUN = function(x) {sum(x==0)/length(x)})
        Variance_estimator = var(Proportion_zero_temps)
        Table_variance_temp = rbind(Table_variance_temp,Variance_estimator)
      }
      mean_estimator = mean(Proportion_zero_temps)
      U = data.frame(Var = log(Table_variance_temp[,1]),
                     N = log(N_cells_values))
      m_var= lm(Var~N,U,subset = is.finite(U$Var))
      predicted_var = exp(predict.lm(object = m_var,newdata = data.frame(N=log(N_cells))))
      return(c(predicted_var,mean_estimator))
    }
    cat(" done ! \n")
  }  
  
  Estimated_variance = approx(x = log10(Mu_values),y = log10(Table_variance[,1]),xout = log10(Mean_expression))
  Estimated_variance = 10^(Estimated_variance$y) # We have the variance of the estimator 
  P_values = pnorm(Proportion_zero,mean = Expected_proportion_zeros,sd = sqrt(Estimated_variance),lower.tail = FALSE)
  P_values_corrected = p.adjust(P_values,method = "fdr")
  
  Selected_genes = names(which(P_values_corrected<P_value_threshold & Delta_excess_zero > Delta_threshold))
  Genes_to_plot = Selected_genes
  pdf(paste0(Output_path,Method,"_",Tissue,"_proportion_zero.pdf"),width = 6.5,height = 6.5,useDingbats = FALSE)
  par(fig=c(0,1,0,1),las=1,bty='l',cex.lab = 0.7, cex.axis = 0.7)
  plot(Total_gene_expression,Proportion_zero,log="x",xlab="Total gene expression",ylab="Proportion 0",
       ylim=c(0,1),yaxs="i",cex.lab=1.15)
  points(Total_gene_expression[order(Total_gene_expression)],Expected_proportion_zeros[order(Total_gene_expression)]
         ,type = "l",lwd=2,col="blue",lty=2)
  points(Total_gene_expression[Selected_genes],Proportion_zero[Selected_genes],pch=21,bg="red")
  text(Total_gene_expression[Genes_to_plot],Proportion_zero[Genes_to_plot],
       labels = Genes_to_plot,pos=3,cex=0.6,offset = 0.1)
  dev.off()
  gc()
  return(list(Selected_genes = Selected_genes,Excess_zero_score = Delta_excess_zero[colnames(Expression_file)]))
}
