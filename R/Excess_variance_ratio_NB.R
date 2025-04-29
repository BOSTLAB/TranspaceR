#' Excess Variance Ratio computation using Negative Binomial Distribution
#'
#' This function analyzes gene expression data to identify genes exhibiting an excess variance ratio.
#' It employs a negative binomial model to calculate p-values and selects significant genes based on
#' user-defined thresholds for p-values and variance ratios and generates a plot of the results..
#'
#' @param Expression_file A data frame or matrix containing gene expression data, where rows represent genes
#' and columns represent samples.
#' @param Method A character string specifying the method used for analysis.
#' @param Tissue A character string indicating the type of tissue being analyzed.
#' @param Output_path A character string specifying the directory where the output PDF will be saved.
#' @param P_value_threshold A numeric value representing the threshold for statistical significance (default is 0.01).
#' @param Ratio_threshold A numeric value representing the threshold for the excess variance ratio (default is 1.5).
#'
#' @return A list containing:
#'      - Selected_genes: A character vector of gene names that show a significant excess variance ratio based on the specified thresholds.
#'      - Excess_variance_ratio: A numeric vector of excess variance ratios for each gene in the expression dataset.
#' @importFrom sparseMatrixStats colVars
#' @import Matrix
#' @import doParallel
#' @import foreach
#' @export

Excess_variance_ratio_NB = function(Expression_file,Output_path,Method,Tissue,P_value_threshold = 0.01,Ratio_threshold = 1.5) {
  
  if (is.data.frame(Expression_file)) {
    Mean_expression = colMeans(Expression_file)
    Var_gene = apply(Expression_file,MARGIN = 2,FUN = var)
  }
  else {
    Mean_expression = Matrix::colMeans(Expression_file)
    Var_gene = sparseMatrixStats::colVars(Expression_file)
  }
  Mean_expression_2 = Mean_expression^2
  m_null= lm((Var_gene) ~ 0 + offset(Mean_expression)+ Mean_expression_2)
  N_cells = nrow(Expression_file)
  N_simulations = 500
  Estimated_s = coef(m_null)
  
  
  Excess_variance_ratio = Var_gene/m_null$fitted.values
  Genes_to_plot = names(Excess_variance_ratio[order(Excess_variance_ratio,decreasing = TRUE)][1:10])
  
  Excess_variance_ratio[is.na(Excess_variance_ratio)] = 0
  
  N_cell_threshold = 10000
  Too_many_cells = FALSE
  
  if (N_cells>N_cell_threshold) {
    Too_many_cells = TRUE
  }
  
  
  Mu_values=quantile(Mean_expression,probs = seq(0,1,length.out = 30))
  if (!Too_many_cells) {
    cat("Computing the regular confidence interval....")
    
    VAriance_estimated <-sapply(Mu_values, function(mu) {
      x = matrix(1, nrow = 1, ncol = N_simulations)
      x = apply(x, MARGIN = 2, FUN = function(x) {
        rnbinom(mu = mu, n = N_cells, size = 1 / Estimated_s)
      })
      Variance_estimated = apply(x, MARGIN = 2, FUN = var)
    })
    
    Table_variance <- sapply(Mu_values, function(mu) {
      x = matrix(1, nrow = 1, ncol = N_simulations)
      x = apply(x, MARGIN = 2, FUN = function(x) {
        rnbinom(mu = mu, n = N_cells, size = 1 / Estimated_s)
      })
      Variance_estimated = apply(x, MARGIN = 2, FUN = var)
      var(Variance_estimated)  # Variance of the variance estimator
    })
    
    rownames(Table_variance) = c()
    cat(" done ! \n")
  }
  
  if (Too_many_cells) {
    cat("Too many cells... confidence interval computed by sub-sampling")
    
    Table_variance <- foreach(i = 1:length(Mu_values), .combine = 'rbind') %dopar% {
      N_cells_values = c(1000, 1500, 2000, 2500, 3000, 5000, 10000)
      
      Variance_estimated <- sapply(N_cells_values, function(n_cells) {
        x = matrix(1, nrow = 1, ncol = N_simulations)
        x = apply(x, MARGIN = 2, FUN = function(x) {
          rnbinom(mu = Mu_values[i], n = n_cells, size = 1 / Estimated_s)
        })
        apply(x, MARGIN = 2, FUN = var)
      })
      
      Table_variance_temp <- sapply(N_cells_values, function(n_cells) {
        x = matrix(1, nrow = 1, ncol = N_simulations)
        x = apply(x, MARGIN = 2, FUN = function(x) {
          rnbinom(mu = Mu_values[i], n = n_cells, size = 1 / Estimated_s)
        })
        Variance_estimated = apply(x, MARGIN = 2, FUN = var)
        var(Variance_estimated)
      })
      mean_estimator = mean(Variance_estimated)
      U = data.frame(Var = log(Table_variance_temp),
                     N = log(N_cells_values))
      m_var = lm(Var ~ N, U, subset = is.finite(U$Var))
      predicted_var = exp(predict.lm(object = m_var, newdata = data.frame(N = log(N_cells))))
      return(c(predicted_var, mean_estimator, summary(m_var)$r.squared))
    }
    cat(" done ! \n")
  }
  
  #We can now say which genes are significantly variable 
  
  Table_variance <- as.data.frame(Table_variance)
  Estimated_variance = approx(x = log10(Mu_values),y = log10(Table_variance[,1]),xout = log10(Mean_expression))
  Estimated_variance = 10^(Estimated_variance$y) # We have the variance of the estimator 
  
  P_values = pnorm(Var_gene,mean = m_null$fitted.values,sd = sqrt(Estimated_variance),lower.tail = FALSE)
  P_values_corrected = p.adjust(P_values,method = "fdr")
  
  Selected_genes = names(which(P_values_corrected<P_value_threshold & Excess_variance_ratio > Ratio_threshold))
  
  Genes_to_plot = Selected_genes
  pdf(paste0(Output_path,Method,"_",Tissue,"_variance.pdf"),width = 6.5,height = 6.5,useDingbats = FALSE)
  par(fig=c(0,1,0,1),las=1,bty='l',cex.lab = 0.7, cex.axis = 0.7)
  plot(Mean_expression,Var_gene,log="xy",xlab="Mean",ylab="Variance",cex.lab=1.3)
  points(Mean_expression[Selected_genes],Var_gene[Selected_genes],pch=21,bg="red")
  points(Mu_values,Mu_values+Mu_values^2*Estimated_s,type='l',lty=2,col="blue",lwd=2)
  text(Mean_expression[Genes_to_plot],Var_gene[Genes_to_plot],
       labels = Genes_to_plot,pos=3,cex=0.6,offset = 0.1)
  title(main = paste0(Method,"_",Tissue), cex.sub = 0.6)
  
  dev.off()
  
  return(list(Selected_genes = Selected_genes, Excess_variance_ratio = Excess_variance_ratio[colnames(Expression_file)]))
}