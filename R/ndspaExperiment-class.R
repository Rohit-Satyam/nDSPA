#' nDSPA
#'
#' Enables easy loading of the data matrices provided by Nanostring.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @importFrom methods setClass new
#' @import dplyr
#' @export
#'
#' @param counts character. Specify the file name in .tsv, .csv or .xlsx
#' @param ... Other arguments passed while creating Summarized/ SingleCell Experiment object.
#' Specify the group file with x,y coordinates (in .csv, .tsv and .xlsx format).
#'
#' @return ne Returns a ndspaExperiment object.

.ndspaEx <- methods::setClass("ndspaExperiment", contains = "SingleCellExperiment")
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Original_ID"))

.ndspaExperiment <- function(counts, ...) {
  ne <- SingleCellExperiment::SingleCellExperiment(list(counts = counts), ...)
  .ndspaEx(ne)
}
