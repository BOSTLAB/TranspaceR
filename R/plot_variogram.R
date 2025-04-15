#' Plot Variogram
#'
#' This function plots the variogram based on the provided parameters and model type.
#'
#' @param Sill Numeric value representing the sill of the variogram.
#' @param Range Numeric value representing the range of the variogram.
#' @param Alpha Numeric value used in the 'Finetuned exponential' model.
#' @param Nugget Numeric value representing the nugget effect.
#' @param Model Character string indicating the model type ('Constant', 'Exponential',
#' 'Finetuned exponential', or 'Cosine exponential').
#' @param filtered_x A data frame containing the variogram data.
#' @param R2 Numeric value representing the coefficient of determination.
#'
#' @return NULL This function does not return a value; it saves a PDF file.
#' @export
#' 
plot_variogram <- function(Sill, Range, Alpha,Gamma, Nugget, Model, filtered_x,R2) {
  par(las=1, bty="l")
  max_x <- max(filtered_x$r) * 1.1
  max_y <- max(filtered_x$Variogram) * 1.2
  plot(filtered_x$r, filtered_x$Variogram, ylim=c(0, max_y), xlim=c(0, max_x),
       xaxs='i', xlab = "r (um)", ylab = "Gamma(r)",
       pch = 19, col = "blue", yaxs='i')
  
  if (Model == 'Constant') {
    abline(h = Sill, col = 'red')
  } else if (Model == 'Exponential') {
    curve(Sill * (1 - exp(-x / Range)) + Nugget, add = TRUE, col = 'red',
          from = 0, to = max_x)
  } else if (Model == 'Finetuned exponential') {
    curve(Sill * (1 - exp(-x^Alpha / Range) ) + Nugget, add = TRUE, col = 'red',
          from = 0, to = max_x)
  }
  legend("topright", legend = paste("R^2 =", round(R2, 3)), col = "red", lwd = 2,
         xjust = 1, yjust = 1, cex = 1)
  return(invisible())
}

