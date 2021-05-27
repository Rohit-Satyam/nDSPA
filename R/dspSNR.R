#' nDSPA
#'
#' This function performs Data Normalization using Signal to Noise Ratio.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom S4Vectors isEmpty
#' @importFrom rlang .data
#' @export
#'
#' @param seobj Provide ndspaExperiment object.
#' @param use Statistical method to be used for scaling.
#' Possible values are "mean" or "gmean".
#' @param probe.set Set of Probes to be used for Normalization.
#' Accepts a vector of custom Probes.
#' Default: "all" uses all probes tagged as "Control" in the rowData.
#' @param use.assay Define on which assay the normalization must be performed.
#' Use assays(seobj) or assayNames(seobj) to view assays saved in ndspaExperiment object.
#' Default:"erccScaled"
#'
#' @return norm_mat numeric. This function returns a normalized matrix.
#' @examples
#' \dontrun{
#' test <- dspSNR(test, use="gmean", probe.set="all", use.assay = "erccScaled")
#'
#' # OR
#'
#' test <- dspSNR(test, use="gmean", probe.set="all",
#' use.assay = "area.scaled.data")
#' }
dspSNR <- function(seobj, use = "gmean", probe.set = "all", use.assay = "erccScaled") {
  probes <- rowData(seobj)

  if (probe.set == "all") {
    Isotype <- probes$`ProbeName (display name)`[probes$`#CodeClass` == "Negative" & probes$`#Analyte type` == "RNA"]
  } else {
    # check if probe.set selected in probe list
    iso_set <- probes$`ProbeName (display name)`[probes$`#CodeClass` == "Negative" & probes$`#Analyte type` == "RNA"]
    not_in_set <- probe.set[!(probe.set %in% iso_set)]
    if (S4Vectors::isEmpty(not_in_set)) {
      cat("All Probes are in Control Probe Set")
    } else {
      warning("Probes not in Isotype probe set: \n", not_in_set, "\n")
    }
    if (all(probe.set %in% iso_set)) {
      Isotype <- probe.set
    } else {
      warning("Not all probes within controls \n")
      Isotype <- setdiff(probe.set, not_in_set)
    }
  }

  Iso_mat <- assay(seobj, use.assay)[Isotype, ]

  if (use == "gmean") {
    if (is.null(dim(Iso_mat))) {
      if (length(Iso_mat) == 0) {
        stop("There are no Control Probes")
      } else if (is.numeric(Iso_mat)) {
        normfactors <- Iso_mat
      } else {
        mode(Iso_mat) <- "numeric"
        normfactors <- Iso_mat
      }
    } else {
      normfactors <- Iso_mat %>%
        apply(2, as.numeric) %>%
        apply(2, log) %>%
        apply(2, mean) %>%
        exp()
    }
  }

  if (use == "mean") {
    if (is.null(dim(Iso_mat))) {
      if (length(Iso_mat) == 0) {
        stop("There are no Control Probes")
      } else if (is.numeric(Iso_mat)) {
        normfactors <- Iso_mat
      } else {
        mode(Iso_mat) <- "numeric"
        normfactors <- Iso_mat
      }
    } else {
      normfactors <- Iso_mat %>%
        apply(2, as.numeric) %>%
        apply(2, mean)
    }
  }
  norm_mat <- t(t(assay(seobj, use.assay)) * (1 / normfactors))
  l <- list(norm_mat)
  names(l) <- paste0(use, "SNR_Normalised")
  assays(seobj) <- c(assays(seobj), l)
  return(seobj)
}


#styler:::style_active_file()
