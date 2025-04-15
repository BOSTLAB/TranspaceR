#' Compute Geary's C Score for Gene Expression Data
#'
#' This function calculates Geary's C score for each gene in the provided expression data,
#' using spatial information from the provided metadata. It also computes p-values and generates a plot of the results.
#'
#' @param Expression_file A matrix or data frame where rows represent samples and columns represent genes.
#' @param Meta_data A data frame containing metadata with `cell_centroid_x` and `cell_centroid_y` representing the spatial coordinates of the cells.
#' @param Method A character string specifying the method used for analysis.
#' @param Tissue A character string indicating the type of tissue being analyzed.
#' @param Output_path A character string specifying the directory where the output PDF will be saved.
#' @param pvalue_threshold A numeric value indicating the threshold for adjusted p-values to select significant genes. Default is 0.01.
#'
#' @return A list containing:
#'    \item{Values}{A named numeric vector of Geary's C scores for each gene.}
#'    \item{Selected_genes}{A character vector of gene names that are significant based on the adjusted p-value threshold.}
#'    \item{P_values}{A numeric vector of adjusted p-values for each gene.}
#' @export
#' 
#' @import igraph
#' @import dplyr
#' @import e1071
Geary_C_score = function(Expression_file,Meta_data,Output_path,Method,Tissue,pvalue_threshold = 0.01 ) {   
  
  temp_delaunay = RCDT::delaunay(cbind(Meta_data$cell_centroid_x,Meta_data$cell_centroid_y))
  graph_delaunay = graph_from_edgelist(el = temp_delaunay$edges[,-3],directed = FALSE )
  spatial_weight_matrix = as_adjacency_matrix(graph_delaunay,sparse = TRUE)
  
  spatial_weight_matrix[is.na(spatial_weight_matrix)] = 0
  spatial_weight_matrix = drop0(spatial_weight_matrix)
  spatial_weight_matrix@x <- 1/spatial_weight_matrix@x #sparse dense coercion
  spatial_weight_matrix = spatial_weight_matrix + t(spatial_weight_matrix)
  
  geary_vector <- c()
  p_value_vector <- c()
  
  S_0 <- sum(spatial_weight_matrix)
  S_1 <- 0.5 * sum((spatial_weight_matrix + t(spatial_weight_matrix))^2)
  x <- colSums(spatial_weight_matrix)
  y <- rowSums(spatial_weight_matrix)
  S_2 <- sum((x + y)^2)
  
  pb = txtProgressBar(min = 1, max = ncol(Expression_file), initial = 1,char = "+",title = "Bla") 
  cat("Computing Geary's C score for each gene \n")
  for (i in 1:ncol(Expression_file)){
    setTxtProgressBar(pb,i)
    X <- matrix(Expression_file[, i], ncol = 1)
    X_squared <- X^2
    
    product_temp_1 <- 2 * sum(spatial_weight_matrix %*% X_squared)
    product_temp_2 <- 2 * t(X) %*% spatial_weight_matrix %*% X
    
    # define parameters for Geary C computation
    N <- length(X)
    W <- sum(spatial_weight_matrix)
    Var_X <- var(X) * N
    
    # compute Geary's C score
    gearys_C <- as.numeric((N - 1) * (product_temp_1 - product_temp_2) / (2 * Var_X * W))
    # append each Geary C score to vector
    geary_vector <- c(geary_vector, gearys_C)
    
    
    b_2 <- kurtosis(Expression_file[,i])
    Var_random <- ((N - 1) * S_1 * (N^2 - 3 * N + 3 - (N - 1) * b_2) -
                     0.25 * (N - 1) * S_2 * (N^2 + 3 * N - 6 - (N^2 - N + 2) * b_2) +
                     S_0^2 * (N^2 - 3 - (N - 1)^2 * b_2)) / (N * (N - 2)^2 * S_0^2)
    p_value = pnorm(gearys_C,mean = 1,sd =sqrt(Var_random),lower.tail = TRUE )
    p_value_vector= c(p_value_vector, p_value)
    
  }
  cat("\n done !")
  
  geary_vector = 1/geary_vector #Geary score = 1/Geary_C
  names(geary_vector) = colnames(Expression_file)
  
  adj_p_value_vector <- p.adjust(p_value_vector,method='fdr')
  names(adj_p_value_vector) = colnames(Expression_file)
  adj_p_value_vector = adj_p_value_vector[order(geary_vector,decreasing = TRUE)]
  
  Values = geary_vector[order(geary_vector,decreasing = TRUE)]
  Selected = names(Values[adj_p_value_vector <= pvalue_threshold])
  
  P_values = adj_p_value_vector
  P_values[P_values<10^-300]=10^-300
  
  pdf(paste0(Output_path,Method,"_",Tissue,"_Geary_score.pdf"),width = 6.5,height = 6.5,useDingbats = FALSE)
  par(las=1, bty='l')
  plot(Values,-log10(P_values),xlab = 'Geary\'s Score', ylab ='P-value (-log10)', main = paste0(Method,"_",Tissue))
  points(Values[Selected],-log10(P_values)[Selected],pch=21,bg="red")
  abline(h=1, lwd=1.5, lty=2, col="red")
  legend('bottomright', legend = c('p-value<0.01', 'p-value>0.01'),
         col = c('red', 'black'), pch = 1, bty = 'l')
  text(Values[1:10],-log10(P_values)[1:10],labels = names(Values[1:10]),pos = 4,cex=0.6)
  dev.off()
  return(list(Values = geary_vector ,Selected_genes = Selected,P_values = P_values))
}

