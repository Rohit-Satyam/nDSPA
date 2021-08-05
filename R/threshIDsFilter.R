#' nDSPA
#'
#' Return IDs to be used in ndspaQC based on various Thresholds.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom rlang .data
#' @export
#'
#' @param seobj S4Vector. Input ndspaExperiment object.
#' @param FOV numeric. Minimum FOV. Default:75 .
#' @param bind.min numeric. Minimum Binding Density. Default: 0.1.
#' @param bind.max numeric. Maximum Binding Density. Default: 2.25.
#' @param min.nuc numeric. Minimum Nucleus Count in AOI. Default:200.
#' @param min.area numeric. Minimum area in AOI. Default: 16000.
#'
#' @return character. Vector of Filtered IDs.
#' @examples
#' \dontrun{
#' ids <- threshIDsFilter(test)
#' }
#'
threshIDsFilter <- function(seobj, FOV = 75, bind.min = 0.1, bind.max = 2.25, min.nuc = 200, min.area = 16000) {
  anno <- data.frame(SummarizedExperiment::colData(seobj), check.names = F)
  out <- anno %>%
    mutate(`Fov counted` = as.numeric(as.character(`Fov counted`))) %>%
    filter("Fov counted" >= FOV) %>%
    mutate(BindingDensity = as.numeric(as.character(BindingDensity))) %>%
    filter(BindingDensity >= bind.min & BindingDensity <= bind.max) %>%
    mutate(`AOI nuclei count` = as.numeric(as.character(`AOI nuclei count`))) %>%
    filter(`AOI nuclei count` >= min.nuc) %>%
    mutate(`AOI surface area` = as.numeric(as.character(`AOI surface area`))) %>%
    filter(`AOI surface area` >= min.area)
  return(out$Original_ID)
}
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Original_ID","AOI surface area",
                                                        "AOI nuclei count","BindingDensity","Fov counted"))
#styler:::style_active_file()
