#' nDSPA
#'
#' This function is DSP_SNR.
#'
#' @import dplyr
#' @export
#'
#' @param df_probe numeric. This is a data frame.
#' @param mat_val_all numeric. This is a data matrix
#' @param method Possible values are "mean" and "geomean". Default Value: "geomean"
#' @param probes Default value is "all"
#'
#' @return norm_mat numeric. This is a data matrix.


DSP_SNR <- function(df_probe,mat_val_all,method="geomean", probes = "all"){
  if (probes == "all") {
    Isotype <- df_probe$`ProbeName (display name)`[df_probe$CodeClass == "Negative" & df_probe$`Analyte type` == "RNA"]
  } else {
    #check if probes selected in probe list
    iso_set <- df_probe$`ProbeName (display name)`[df_probe$CodeClass == "Negative" & df_probe$`Analyte type` == "RNA"]
    not_in_set <- probes[!(probes %in% iso_set)]
    cat("Probes not in Isotype probe set: \n")
    print(not_in_set)
    if (all(probes %in% iso_set)) {
      Isotype <- probes
    } else {
      cat("Not all probes within controls \n")
      Isotype <- probes
    }
  }
  Iso_mat <- mat_val_all[Isotype,]

  if (method=="geomean") {
    if (is.null(dim(Iso_mat))){
      if(length(Iso_mat)==0){
        print("There are no Control Probes") #Change to error condition
      }else if (is.numeric(Iso_mat)){
        normfactors <- Iso_mat
      }else{
        mode(Iso_mat) <- "numeric"
        normfactors <- Iso_mat
      }
    }else{
      normfactors <- Iso_mat %>%
        apply(2, as.numeric) %>%
        apply(2, log) %>%
        apply(2, mean) %>%
        exp()
    }


  }

  if (method=="mean") {
    if (is.null(dim(Iso_mat))){
      if(length(Iso_mat)==0){
        print("There are no Control Probes") #Change to error condition
      }else if (is.numeric(Iso_mat)){
        normfactors <- Iso_mat
      }else{
        mode(Iso_mat) <- "numeric"
        normfactors <- Iso_mat
      }
    }else{
      normfactors <- Iso_mat %>%
        apply(2, as.numeric) %>%
        apply(2, mean)
    }


  }
  norm_mat <- t(t(mat_val_all)*(1/normfactors))

  return(norm_mat)
}

