#' nDSPA
#'
#' This function performs Scaling.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom rlang .data
#' @export
#'
#' @param seobj Provide ndspaExperiment object.
#' @param use Statistical method to be used for scaling. Possible values are "mean" or "gmean".
#' @param method Variable used to scale data. Possible values are "area" and "nuclei". Default: "area".
#'
#' @return seobj Returns ndspaExperiment updated with scaled assay.
#' @examples
#' \dontrun{
#' test <- ndspaScale(test,use ="gmean", method = "area")
#'
#' # OR
#'
#' test <- ndspaScale(test,use ="gmean", method = "nuclei")
#' }

ndspaScale <- function(seobj, use = c("mean", "gmean"), method = "area") {
  anno <- SummarizedExperiment::colData(seobj)
  counts <- SummarizedExperiment::assay(seobj, "erccScaled")
  if (!tolower(method) %in% c("area", "nuclei")) {
    stop("Error in method: Must be area or nuclei\n")
  }

  if (use == "gmean") {
    if (method == "area") {
      geomean_scale_area <- (anno$`AOI surface area` %>% as.character() %>% as.numeric() %>% gmean()) / (anno$`AOI surface area` %>% as.character() %>% as.numeric())
      scaled_df <- t(t(counts) * geomean_scale_area)
    } else {
      geomean_scale_nuclei <- (anno$`AOI nuclei count` %>% as.character() %>% as.numeric() %>% gmean()) / (anno$`AOI nuclei count` %>% as.character() %>% as.numeric())
      scaled_df <- t(t(counts) * geomean_scale_nuclei)
    }
    l <- list(scaled_df)
    names(l) <- "scaled.data"
    assays(seobj) <- c(assays(seobj), l)
    return(seobj)
  } else if (use == "mean") {
    if (method == "area") {
      mean_scale_area <- (anno$`AOI surface area` %>% as.character() %>% as.numeric() %>% mean()) / (anno$`AOI surface area` %>% as.character() %>% as.numeric())
      scaled_df <- t(t(counts) * mean_scale_area)
    } else {
      mean_scale_nuclei <- (anno$`AOI nuclei count` %>% as.character() %>% as.numeric() %>% mean()) / (anno$`AOI nuclei count` %>% as.character() %>% as.numeric())
      scaled_df <- t(t(counts) * mean_scale_nuclei)
    }
    l <- list(scaled_df)
    names(l) <- paste0(method,".scaled.data")
    assays(seobj) <- c(assays(seobj), l)
    return(seobj)
  } else {
    stop("Provide use argument.")
  }
}
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "gmean", "y"))
#styler:::style_active_file()
