#' nDSPA
#'
#' This function performs Scaling steps
#'
#' @import dplyr
#' @export
#'
#' @param df_annonumeric. This is a data frame
#' @param mat_val_all numeric. This is a data matrix.
#' @param method Possible values are "area" and "nuclei"
#'
#' @return scaled_df numeric. This is a dataframe.


DSP_Scale <- function(df_anno, mat_val_all,method="area"){
  if (!tolower(method) %in% c("area","nuclei")){
    cat("Error in method: Must be area or nuclei\n")
  } else if(method=="area") {
    geomean_scale_area <- (df_anno$`AOI surface area` %>% as.character() %>% as.numeric() %>% gmean())/(df_anno$`AOI surface area` %>% as.character() %>% as.numeric())
    scaled_df <- t(t(mat_val_all)*geomean_scale_area)
  } else {
    geomean_scale_nuclei <- (df_anno$`AOI nuclei count` %>% as.character() %>% as.numeric() %>% gmean())/(df_anno$`AOI nuclei count` %>% as.character() %>% as.numeric())
    scaled_df <- t(t(mat_val_all)*geomean_scale_nuclei)
  }
  return(scaled_df)
}
