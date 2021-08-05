.gmean <- function(x,method="log") {
  if (method == "log"){
    #Safer method does not produce overflows
    gm <- exp(mean(log(x)))

  }else if (method == "mult"){
    gm <- prod(x)^(1/len(x))
  }
  return(gm)
}
