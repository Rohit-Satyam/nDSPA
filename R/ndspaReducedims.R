#' nDSPA
#'
#' Adding low-dimensional representations to nDSPA object
#'
#' The ndspaReducedims function let user perform Dimension Reduction given ndspaExperiment
#' object and returns a matrix with computed PCs or other low dimensional embedding depending on type of method used.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import M3C
#' @importFrom stats prcomp
#' @import plotly
#' @import ggplot2
#' @export
#'
#' @param seobj S4Vector Provide ndspaExperiment object.
#' @param use.assay character Provide the assay name to be used for QC plots. Default "counts".
#' @param type character Type of Dimensional Reduction method to be used ("PCA","tsne","umap"). Default: "PCA"
#' @param ... arguments passed to or from prcomp or Rtsne depending on type of method opted for dimension reduction.
#' @return Matrix containing the new representations for the ndspaExperiment object.
#' @examples
#' \dontrun{
#'
#' # Calculating PCA
#' test.pca <- ndspaReducedims(test)
#'
#' # Calculating tSNE
#' test.tsne <- ndspaReducedims(test, type = "tsne")
#'
#' #Calculating UMAP
#'
#' test.umap <- ndspaReducedims(test, type = "umap")
#' # Saving results
#' reducedDims(test) <- list(PCA = test.pca$data, TSNE = test.tsne$data, UMAP=test.umap$data)
#' }
#'
ndspaReducedims <- function(seobj, use.assay = "counts", type = "PCA", ...) {
  anno <- data.frame(SummarizedExperiment::colData(seobj), check.names = FALSE)
  val_all <- assay(seobj, assay = use.assay)
  probe <- data.frame(rowData(seobj), check.names = FALSE)
  Endo_probes <- probe$`ProbeName (display name)`[probe$`#CodeClass` == "Endogenous"]
  val_Endo <- val_all[rownames(val_all) %in% Endo_probes, ]%>% apply(2, log)

  if (type == "tsne") {
    # tsne_data <- Rtsne::Rtsne(pca_data$x[, seq_len(rank)], perplexity = perplexity, pca = FALSE, ...)
    tsne_data <- M3C::tsne(val_Endo,dotsize = 3,...)+ggplot2::theme_minimal()
    rownames(tsne_data$data) <- colnames(val_Endo)
    return(tsne_data)
  } else if (type == "PCA") {
    # pca_data <- M3C::umap(val_Endo,dotsize = 3,...)+ggplot2::theme_minimal()
    pca_data <- stats::prcomp(t(val_Endo), ...)
    return(pca_data)
  } else if (type == "umap"){
    umap_data <- M3C::umap(val_Endo,dotsize = 3,...)+ggplot2::theme_minimal()
    return(umap_data)
  }else {
    stop("Enter appropriate dimension reduction method using argument: type")
  }
}

if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Scan_ID", "#CodeClass", "#Analyte type", "ID", "ProbeName (display name)", "Segment tags", "ROI_ID", "group", "len", "row.name"))
# styler:::style_active_file()
