#' Convert Strings to Colors
#'
#' This function maps unique strings to colors. In absence of a provided color palette,
#' it will generate a rainbow palette.
#'
#' @param string A character vector or factor containing the strings to be converted
#'                            to colors.
#' @param colors A character vector of colors to be used for the unique elements
#'                            in `string`. If NULL, a rainbow palette will be generated.
#' @return A character vector of colors corresponding to the input strings.
#' 
#' @examples
#' annotation = c('B cells', 'T cells' , 'B cells', 'pDCs')
#' colors = string.to.colors(annotation)
#' print(colors)
#' [1] "#FF0000" "#00FF00" "#FF0000" "#0000FF"
#' 
#' @export

string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}

