#' nDSPA
#'
#' This function performs basic mean.
#'
#' @import dplyr
#' @export
#'
#' @param x numeric. This is a data matrix.
#' @param method numeric. The possible values are "log" and "mult". Default Value: "log"
#'
#' @return gm numeric. Returns a gmeans


to.numeric <- function(x) as.numeric(as.character(x))
gmean <- function(x,method="log") {
  if (method == "log"){
    #Safer method does not produce overflows
    gm <- exp(mean(log(x)))

  }else if (method == "mult"){
    gm <- prod(x)^(1/len(x))
  }
  return(gm)
}
