#' nDSPA
#'
#' Return IDs to be used in ndspaQC.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom rlang .data
#' @export
#'
#' @param seobj S4Vector. Input ndspaExperiment object.
#' @param  PCF.min numeric. Provide a numeric threshold.
#' @param PCF.max numeric. Provide a numeric threshold.
#' @return filtered_scalefactor Vector of Filtered IDs.
#' @examples
#' \dontrun{
#' ids <- pcfIDsFilter(test)
#' }

pcfIDsFilter <- function(seobj, PCF.min = 0.3, PCF.max = 3.0) {
  anno <- SummarizedExperiment::colData(seobj)
  scalefactor <- SummarizedExperiment::colData(seobj)$scalefactor
  names(scalefactor) <- SummarizedExperiment::colData(seobj)$Original_ID
  filtered_scalefactor <- scalefactor[scalefactor >= PCF.min & scalefactor <= PCF.max]
  filtred_out_by_posctrl <- names(scalefactor[scalefactor <= PCF.min | scalefactor >= PCF.max])
  ## Or use this if (length(names(scalefactor))==length(filtered_scalefactor))
  if (all(names(scalefactor)%in%names(filtered_scalefactor))) {
    warning("No Samples filtred out by Positive Control Factor \n")
  } else {
    warning(paste0("Failed Positive Control Factor: \n", paste(filtred_out_by_posctrl, collapse = ", "), "\n"))
  }
  return(names(filtered_scalefactor))
}

#styler:::style_active_file()
