#' nDSPA
#'
#' This function performs Quality Control steps
#'
#' @import dplyr
#' @export
#'
#' @param anno numeric. This is a dataframe.
#' @param val_all_df numeric. This is a dataframe
#' @param scalefactor numeric.
#' @param thresh_filt logical. Possible values are TRUE or FALSE.
#' @param PCF_filt logical. Possible values are TRUE or FALSE.


DSP_QC <- function(anno, val_all_df, scalefactor, thresh_filt=FALSE, PCF_filt=FALSE){
  if (!isFALSE(thresh_filt)){
    ID <- NULL
    anno <- anno %>% filter(ID %in% thresh_filt)
  }
  if (!isFALSE(PCF_filt)){
    anno <- anno %>% filter(ID %in% PCF_filt)
  }

  sf <- scalefactor[anno$ID]
  vadf <- val_all_df[,anno$ID]
  if (all.equal(names(scalefactor), colnames(val_all_df))){
    QCdf <- t(t(vadf)*sf)
  }else{
    cat("Scale values not in filtered data: /n")
    print(names(scalefactor[!names(scalefactor) %in% anno$ID]))

    QCdf <- t(t(vadf)*sf)
  }
  return(QCdf) ##<< this is the DF not the anno
}
