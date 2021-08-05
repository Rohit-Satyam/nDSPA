#' nDSPA
#'
#' This function performs Data Normalization steps
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom S4Vectors isEmpty
#' @export
#'
#' @param seobj Provide ndspaExperiment object.
#' @param use Statistical method to be used for scaling.
#' Possible values are "mean" or "gmean".
#' @param probe.set Set of Probes to be used for Normalization.
#' Accepts a vector of custom Probes.
#' Default: "all" uses all probes tagged as "Control" in the rowData.
#' @param use.assay Define on which assay the normalization must be performed.
#' Use assays(seobj) or assayNames(seobj) to view assays saved in ndspaExperiment
#' object. Default:"erccScaled"
#'
#' @return norm_mat numeric. This function returns a normalized matrix.
#' @examples
#' \dontrun{
#' test <- dspNormalization(test, use = "gmean", probe.set="all",
#' use.assay = "erccScaled")
#' # OR
#' test <- dspNormalization(test, use = "gmean", probe.set="all",
#' use.assay = "areaScaled")
#' }

dspNormalization <- function(seobj, use = "gmean", probe.set = "all", use.assay = "erccScaled") {
  probes <- SummarizedExperiment::rowData(seobj)
  if (isTRUE(probe.set == "all")) {
    Controls <- probes$`ProbeName (display name)`[probes$`#CodeClass` == "Control" & (probes$`#Analyte type` == "RNA" | probes$`#Analyte type` == "Protein")]
  } else {
    # check if probes selected in probe list
    ctrl_set <- probes$`ProbeName (display name)`[probes$`#CodeClass` == "Control" & (probes$`#Analyte type` == "RNA" | probes$`#Analyte type` == "Protein")]
    not_in_set <- probe.set[!(probe.set %in% ctrl_set)]
    if (S4Vectors::isEmpty(not_in_set)) {
      cat("All Probes are in Control Probe Set")
    } else {
      warning("Probes not in control probe set: \n", not_in_set, "\n")
    }
    if (all(probe.set %in% ctrl_set)) {
      Controls <- probe.set
    } else {
      warning("Not all probes within controls \n")
      Controls <- setdiff(probe.set, not_in_set)
    }
  }
  Ctrl_mat <- assay(seobj, use.assay)[Controls, ]

  if (use == "gmean") {
    if (is.null(dim(Ctrl_mat))) {
      if (length(Ctrl_mat) == 0) {
        stop("There are no Control Probes")
      } else if (is.numeric(Ctrl_mat)) {
        normfactors <- Ctrl_mat
      } else {
        mode(Ctrl_mat) <- "numeric"
        normfactors <- Ctrl_mat
      }
    } else {
      normfactors <- Ctrl_mat %>%
        apply(2, as.numeric) %>%
        apply(2, log) %>%
        apply(2, mean) %>%
        exp()
    }

    scalefactor <- mean(normfactors) / normfactors
  }

  if (use == "mean") {
    if (is.null(dim(Ctrl_mat))) {
      if (length(Ctrl_mat) == 0) {
        stop("There are no Control Probes")
      } else if (is.numeric(Ctrl_mat)) {
        normfactors <- Ctrl_mat
      } else {
        mode(Ctrl_mat) <- "numeric"
        normfactors <- Ctrl_mat
      }
    } else {
      normfactors <- Ctrl_mat %>%
        apply(2, as.numeric) %>%
        apply(2, mean)
    }

    scalefactor <- mean(normfactors) / normfactors
  }
  norm_mat <- t(t(assay(seobj, use.assay)) * scalefactor)
  l <- list(norm_mat)
  names(l) <- paste0(use, "HK_Normalised")
  assays(seobj) <- c(assays(seobj), l)
  S4Vectors::metadata(seobj) <- c(S4Vectors::metadata(seobj), "Normalization" = paste("Housekeeping normalization was performed on",use.assay,"using",use,". Probset used were:",probe.set))
  return(seobj)
}

#styler:::style_active_file()
