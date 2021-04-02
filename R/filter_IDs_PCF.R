#' nDSPA
#'
#' This function filter IDs by positive control Scaling
#'
#' @import dplyr
#' @export
#'
#' @param anno numeric. This is a dataframe.
#' @param scalefactor numeric. This is data frame.
#' @param PCF_min numeric. Default Value: 0.3
#' @param PCF_max numeric. Default value: 3.0
#'
#' @return filtered_scalefactor numeric. This is a data frame.


#========nDSPA ====
#returns names of IDs which pass filter check
filter_IDs_PCF <- function(anno, scalefactor, PCF_min=0.3, PCF_max=3.0){
  filtered_scalefactor <- scalefactor[scalefactor >= PCF_min & scalefactor<= PCF_max]

  filtred_out_by_posctrl <- names(scalefactor[scalefactor <= PCF_min & scalefactor>= PCF_max])
  if(length(filtred_out_by_posctrl) == 0){
    cat("No Samples filtred out by Positive Control Factor \n")
  }else{
    cat(paste0("Failed Positive Control Factor: \n", paste(names(scalefactor[scalefactor >= PCF_min & scalefactor<= PCF_max]),collapse = ", "),"\n"))
  }
  return(names(filtered_scalefactor))
}
