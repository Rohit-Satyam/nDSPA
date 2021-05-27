#' nDSPA
#'
#' This function performs Quality Control steps.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom S4Vectors metadata
#' @importFrom rlang .data
#' @export
#'
#' @param seobj S4Vector. Input ndspaExperiment object.
#' @param thresh_filt vector. Provide Segment displayed name IDs as a vector.
#' @param PCF_filt vector. Provide Segment displayed name IDs as a vector.
#' @param tag.only logical. To tag the cells that passes default QC criteria. Possible values are TRUE or FALSE.
#' @return seobj Updated ndspaExperiment object with erccScaled assay added.
#' @examples
#' \dontrun{
#' # Run without filtering
#' test <- ndspaQC(test)
#'
#' # Run with PCF Filter ids
#' test <- ndspaQC(test,PCF_filt=ids)
#'
#' # Run in Tag only mode
#' test <- ndspaQC(test,PCF_filt=ids,tag.only=TRUE)
#'
#' }

ndspaQC <- function(seobj, thresh_filt = NULL, PCF_filt = NULL, tag.only = FALSE) {
  if (isFALSE(tag.only)) {
    anno <- SummarizedExperiment::colData(seobj)
    if (!is.null(thresh_filt)) {
      anno <- data.frame(anno, check.names = FALSE)
      anno <- anno %>% dplyr::filter(.data$Original_ID %in% thresh_filt)
    }
    if (!is.null(PCF_filt)) {
      anno <- data.frame(anno, check.names = FALSE)
      anno <- anno %>% dplyr::filter(Original_ID %in% PCF_filt)
      S4Vectors::metadata(seobj) <- c(S4Vectors::metadata(seobj), "QC.Filter.IDs" = list(PCF_filt))
    }

    if (isFALSE("scalefactor" %in% colnames(colData(seobj)))) {
      stop("ScaleFactors are missing. Run erccScaleFactor() function first")
    }
    ## recreating scalefactor obj:########
    scalefactor <- SummarizedExperiment::colData(seobj)$scalefactor
    names(scalefactor) <- colnames(seobj)
    ######################################
    sf <- scalefactor[anno$Original_ID]

    vadf <- SummarizedExperiment::assay(seobj)[, anno$Original_ID]
    if (all.equal(names(scalefactor), colnames(seobj))) {
      QCdf <- t(t(vadf) * sf)
    } else {
      warning("Scale values not in filtered data: /n")
      print(names(scalefactor[!names(scalefactor) %in% anno$Original_ID]))

      QCdf <- t(t(vadf) * sf)
    }
    seobj <- seobj[, colnames(QCdf)]
    l <- list(QCdf)
    names(l) <- "erccScaled"
    SummarizedExperiment::assays(seobj) <- c(SummarizedExperiment::assays(seobj), l)
    return(seobj) ## << this is the DF not the anno
  } else if (!isFALSE(tag.only) & is.null(PCF_filt)) {
    stop("No probe IDs provided using PCF_filt option.")
  } else {
    temp <- gsub("TRUE", "Retained", rownames(SummarizedExperiment::colData(seobj)) %in% PCF_filt)
    temp <- gsub("FALSE", "Removed", temp)
    SummarizedExperiment::colData(seobj)$tag <- temp
  }
  return(seobj)
}

if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Original_ID"))
#styler:::style_active_file()
