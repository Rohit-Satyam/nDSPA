#' nDSPA
#'
#' This function performs Data Normalization steps
#'
#' @import dplyr
#' @export
#'
#' @param df_probe numeric
#' @param mat_val_all numeric
#' @param method Possible methods of normalization include "geomean" and "mean".
#' @param probes The probes argument checks the presence of all probes within controls.
#'
#' @return norm_mat numeric. This function returns a normalized matrix.


DSP_normalization <- function(df_probe,mat_val_all,method="geomean", probes="all"){
  if (probes == "all") {
    Controls <- df_probe$`ProbeName (display name)`[df_probe$CodeClass == "Control" & df_probe$`Analyte type` == "RNA"]
  } else {
    #check if probes selected in probe list
    ctrl_set <- df_probe$`ProbeName (display name)`[df_probe$CodeClass == "Control" & df_probe$`Analyte type` == "RNA"]
    not_in_set <- probes[!(probes %in% ctrl_set)]
    cat("Probes not in control probe set: \n")
    print(not_in_set)
    if (all(probes %in% ctrl_set)) {
      Controls <- probes
    } else {
      cat("Not all probes within controls \n")
      Controls <- probes
    }
  }
  Ctrl_mat <- mat_val_all[Controls,]

  if (method=="geomean") {
    if (is.null(dim(Ctrl_mat))){
      if(length(Ctrl_mat)==0){
        print("There are no Control Probes") #Change to error condition
      }else if (is.numeric(Ctrl_mat)){
        normfactors <- Ctrl_mat
      }else{
        mode(Ctrl_mat) <- "numeric"
        normfactors <- Ctrl_mat
      }
    }else{
      normfactors <- Ctrl_mat %>%
        apply(2, as.numeric) %>%
        apply(2, log) %>%
        apply(2, mean) %>%
        exp()
    }

    scalefactor <- mean(normfactors)/normfactors
  }

  if (method=="mean") {
    if (is.null(dim(Ctrl_mat))){
      if(length(Ctrl_mat)==0){
        print("There are no Control Probes") #Change to error condition
      }else if (is.numeric(Ctrl_mat)){
        normfactors <- Ctrl_mat
      }else{
        mode(Ctrl_mat) <- "numeric"
        normfactors <- Ctrl_mat
      }
    }else{
      normfactors <- Ctrl_mat %>%
        apply(2, as.numeric) %>%
        apply(2, mean)
    }

    scalefactor <- mean(normfactors)/normfactors
  }
  norm_mat <- t(t(mat_val_all)*scalefactor)

  return(norm_mat)
}
