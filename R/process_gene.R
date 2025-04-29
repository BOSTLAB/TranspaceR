#' Process to compute the Variogram of one Gene
#'
#' This function computes variogram for a gene and tries fitting different models to the variogram and 
#' selects the best model based on the Akaike Information Criterion (AIC).
#'
#' @param gene_name A string representing the name of the gene to be processed.
#' @param sum_by_cut A data frame containing the rasterized expression data.
#' @param n_pad A numeric value used to determine the size of unique elements in the rasterized data.
#' @param save_plot A logical value indicating whether to save the plot of the variogram.
#' @param Output_path A string specifying the path where output files should be saved.
#'
#' @return A numeric vector containing the following elements:
#' \item{Sill}{The sill of the best fitting model.}
#' \item{Range}{The range of the best fitting model (if applicable).}
#' \item{Nugget}{The nugget effect of the best fitting model (if applicable).}
#' \item{Alpha}{The alpha parameter of the best fitting model (if applicable).}
#' \item{Total_variance}{The total variance of the best fitting model.}
#' \item{R2}{The R-squared value of the best fitting model.}
#' \item{Prop}{The proportion of variance explained by the spatial structure of the best fitting model.}
#' \item{AIC}{The AIC value of the best fitting model.}
#' \item{Model}{A string indicating the name of the best fitting model.}
#' \item{Variogram}{A list of variogram values used in the fitting process.}
#' @export
#' @importFrom reshape2 dcast
#' @import stats
#' @import dplyr
#' 
process_gene <- function(gene_name,sum_by_cut,n_pad,save_plot,Output_path) {
  #matrix creation
  xl <- reshape2::dcast(sum_by_cut, portion_x ~ portion_y, value.var = gene_name, fun.aggregate = mean)
  xl <- as.matrix(xl[,-1])
  image(xl)
  x = Get_isotropic_vario(Variogram=Get_variogram_map(xl))
  x = x[x$r > 0, ]
  filtered_x <- x[x$r < 30, ]
  
  filtered_x <- filtered_x %>%
    mutate(group = (row_number() - 1) %/% 10) %>%
    group_by(group) %>%
    summarise(
      Variogram = mean(Variogram, na.rm = TRUE),
  r = max(r, na.rm = TRUE))
  filtered_x$r = filtered_x$r * n_pad
  #fitting
  n <- min(unlist(filtered_x$Variogram))
  tau <- filtered_x$r[which.max(filtered_x$Variogram)]
  alpha = 1
  a = filtered_x$r
  b = filtered_x$Variogram
  
  AIC0 <- AIC1 <- AIC2 <- NA
  m_0 <- m_1 <- m_2  <- NA
  Sill<-Range<-Nugget<-Alpha <-Total_variance<-R2<-Prop<-AIC<-NA
  m_0 <- tryCatch({
    stats::nls(b ~ C, start = list(C = mean(b, na.rm = TRUE)))
  }, error = function(e) {
    #  message("Fitting by constant model failed")
  })
  
  if (!all(is.na(m_0))) {
    predicted_values0 <- predict(m_0)
    error2 <- sum((b-predicted_values0)^2)
    AIC0 <- length(b) * log(error2/length(b)) + 2 * 1
    
  }
  
  m_1 <- tryCatch({stats::nls(formula = b ~ C * (1 - exp(-a/tau)) + n,
                              start = list(tau = tau, n = n, C = max(filtered_x$Variogram)-n),lower = c(0,0,0),algorithm = "port")
  } , error = function(e) {
    #message("Fitting by exponential model failed")
  })
  
  if (!all(is.na(m_1))) {
    predicted_values1 <- predict(m_1)
    error2 <- sum((b-predicted_values1)^2)
    AIC1 <- length(b) * log(error2/length(b)) + 2 *3 
  }
  
  m_2 <- tryCatch({stats::nls(formula = b ~ C * (1 - exp(-a^alpha/tau)) + n,
                              start = list(tau = tau, n = n, C = max(filtered_x$Variogram)-n,alpha=alpha),lower = c(0,0,0,0),algorithm = "port")
  } , error = function(e) {
    #message("Fitting by finetuned exponential model failed")
  })
  if (!all(is.na(m_2))) {
    predicted_values2 <- predict(m_2)
    error2 <- sum((b-predicted_values2)^2)
    AIC2 <- length(b) * log(error2/length(b)) + 2 *4 
  }
  
  # Determine the model with the smallest AIC
  AIC_values <- c(AIC0, AIC1, AIC2)
  best_model_index <- which.min(AIC_values)
  AIC = AIC_values[best_model_index]
  
  if (best_model_index == 1) {
    best_model <- m_0
    Sill <- coef(best_model)[1]  
    Model <- 'Constant'
    R2 <- 1 - (sum((filtered_x$Variogram - predict(best_model))^2) / sum((filtered_x$Variogram - min(filtered_x$Variogram, na.rm = TRUE))^2))
  } else if (best_model_index == 2) {
    best_model <- m_1
    Sill <- coef(best_model)[3]
    Range <- coef(best_model)[1]
    Nugget <- coef(best_model)[2]
    Total_variance <- coef(best_model)[3] + coef(best_model)[2]
    Prop <- coef(best_model)[3] / (coef(best_model)[2] + coef(best_model)[3])
    Model <- 'Exponential'
    R2 <- 1 - (sum((filtered_x$Variogram - predict(best_model))^2) / sum((filtered_x$Variogram - mean(filtered_x$Variogram, na.rm = TRUE))^2))
  } else if (best_model_index == 3) {
    best_model <- m_2
    Sill <- coef(best_model)[3]
    Nugget <- coef(best_model)[2]
    Alpha <- coef(best_model)[4]
    Range <- coef(best_model)[1]
    Total_variance <- coef(best_model)[3] + coef(best_model)[2]
    Prop <- coef(best_model)[3] / (coef(best_model)[2] + coef(best_model)[3])
    Model <- 'Finetuned exponential'
    R2 <- 1 - (sum((filtered_x$Variogram - predict(best_model))^2) / sum((filtered_x$Variogram - mean(filtered_x$Variogram, na.rm = TRUE))^2))
   } 
  if (save_plot == TRUE){
    pdf(paste0(Output_path,"Variogram_plots/",gene_name,".pdf"), useDingbats = FALSE, width = 8, height = 8)
    plot_variogram(Sill,Range,Alpha,Gamma,Nugget,Model,filtered_x,R2)
    dev.off()
  }
  return(list(Sill = Sill, Range = Range, Nugget = Nugget, Alpha = Alpha,
              Total_variance = Total_variance, R2 = R2,
              Prop = Prop, AIC = AIC, Model = Model, Vario_values = c(filtered_x$Variogram)))
}

