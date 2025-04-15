#' Compute the radius of cells based on the specified method.
#'
#' @param Meta_data A data frame containing cell metadata including area or volume.
#' @param Method A string specifying the method for radius calculation ('Xenium', 'CosMx', or 'Merfish').
#' @return The updated Meta_data with calculated radius values.
#' @export
#' 
Compute_radius = function(Meta_data, Method) {
  if (Method=='Xenium'){
    Meta_data$radius = sqrt(Meta_data$Area/pi)  
  } else if (Method=='CosMx'){
    Meta_data$radius = sqrt((Meta_data$Area/100)/pi)
  } else if  (Method=='Merfish')  {
    Meta_data$radius = ((3 * Meta_data$volume) / (4 * pi))**(1/3)
  } else {
    'Please choose one method between Xenium, CosMx and Merfish'
  }
  return(Meta_data)
}
