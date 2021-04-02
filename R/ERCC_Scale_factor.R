#' nDSPA
#'
#' This function performs data scaling steps
#'
#' @import dplyr
#' @export
#'
#' @param df_probe numeric This is input.
#' @param df_val_all numeric Input the datamatrix.
#' @return scalefactor numeric.
#'


ERCC_Scale_factor <- function(df_probe, df_val_all){
  ERCC_Probes <- df_probe$`ProbeName (display name)`[df_probe$CodeClass == "Positive" & df_probe$`Analyte type` == "SpikeIn"]

  #check before here that input values are data.matrix with numeric data
  PosCtrl_mat <- df_val_all[ERCC_Probes,]

  #Need try catch handling
  if (is.null(dim(PosCtrl_mat))){
    if(length(PosCtrl_mat)==0){
      print("There are no ERCC_Probes") #Change to error condition
    }else if (is.numeric(PosCtrl_mat)){
      normfactors <- PosCtrl_mat
    }else{
      mode(PosCtrl_mat) <- "numeric"
      normfactors <- PosCtrl_mat
    }
  }else{
    normfactors <- PosCtrl_mat %>%
      apply(2, as.numeric) %>%
      apply(2, log) %>%
      apply(2, mean) %>%
      exp()
  }

  scalefactor <- mean(normfactors)/normfactors

  return(scalefactor)
}
