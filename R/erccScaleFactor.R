#' nDSPA
#'
#' Enables calculation of ERCC ScaleFactors.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @import dplyr
#' @importFrom rlang .data
#' @export
#'
#' @param seobj S4Vector Provide ndspaExperiment object.
#' @return seobj ndspaExperiment Object with ColData populated with ERCC ScaleFactors.
#' @examples
#' \dontrun{
#' test <- erccScaleFactor(test)
#' }

erccScaleFactor <- function(seobj) {
  probes <- SummarizedExperiment::rowData(seobj) # DFrame
  val.all <- SummarizedExperiment::assay(seobj) ## Matrix object
  ERCC_Probes <- probes$`ProbeName (display name)`[probes$`#CodeClass` == "Positive" & probes$`#Analyte type` == "SpikeIn"]
  # check before here that input values are data.matrix with numeric data
  PosCtrl_mat <- val.all[ERCC_Probes, ]


  # trycatch handling
  if (is.null(dim(PosCtrl_mat))) {
    if (length(PosCtrl_mat) == 0) {
      stop("There are no ERCC_Probes")
    } else if (is.numeric(PosCtrl_mat)) {
      normfactors <- PosCtrl_mat
    } else {
      mode(PosCtrl_mat) <- "numeric"
      normfactors <- PosCtrl_mat
    }
  } else {
    normfactors <- PosCtrl_mat %>%
      apply(2, as.numeric) %>%
      apply(2, log) %>%
      apply(2, mean) %>%
      exp()
  }

  scalefactor <- mean(normfactors) / normfactors
  SummarizedExperiment::colData(seobj)$scalefactor <- data.frame(scalefactor)$scalefactor
  S4Vectors::metadata(seobj) <- c(S4Vectors::metadata(seobj), "erccScaleFactor" = "ERCC Scale Factors were computed.")
  return(seobj)
}

#styler:::style_active_file()
