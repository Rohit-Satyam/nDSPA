#' nDSPA
#'
#' Enables easy loading of the data matrices provided by Nanostring.
#'
#' The readnDSPA function is built at top of SingleCellExperiment object and is used to load Nanostring spatial omics data.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom S4Vectors merge metadata
#' @importFrom rlang .data
#' @importFrom readr read_tsv
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom tools file_ext
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @param x character. Specify the file name in .tsv, .csv or .xlsx
#' @param experiment character. Specify the experiment type whether it is "Protein", "RNA" or "CTA".
#' @param meta This is an experimental argument.
#' Specify the group file with x,y coordinates (in .csv, .tsv and .xlsx format).
#'
#' @return ne Returns a ndspaExperiment object.
#' @examples
#' \dontrun{
#' test <- readnDSPA("01-1.dsp_data.raw.sim.tsv", experiment = "RNA",
#' meta = "01-2.dsp_roi.metadata.sim.tsv")
#' }

readnDSPA <- function(x, experiment = c("Protein", "RNA", "CTA"), meta = NULL) {
  ####### Importing Data ####################
  inputext <- tools::file_ext(x)
  if (!is.null(meta)) {
    metaext <- tools::file_ext(meta)
    ifelse(metaext == "tsv" | metaext == "csv" | metaext == "xlsx",
           print("Reading Metadata"),
           stop("The selected file type for metadata is not available.
                Metadata must be in csv,tsv or xlsx format."))
    meta_df <- .readfn(metaext, meta) %>%
      dplyr::mutate(names = paste(.$Scan_ID, .$ROI, .$segment, sep = " | ")) %>%
      `rownames<-`(.$names)
  }
  # Imp function/test data type

  df <- .readfn(inputext, x)



  ################# Experiment type#####################################
  experiment <- experiment

  ifelse(experiment == "Protein" | experiment == "RNA" | experiment ==
           "CTA", exp_type <- list(experiment = "DataType"),
         stop("Experiment type not provided"))
  names(exp_type) <- experiment
  ###### Parsing Data and making ndspa Object##########################

  anno <- .cut_anno(df)
  probes <- .cut_probes(df)
  val.All <- .cut_vals(df, anno, probes)

  if (!is.null(meta)) {
    anno <- S4Vectors::merge(x = anno, y = select(meta_df, c(names,x, y)),
                             by.x = "Original_ID", by.y = "names", all.x = TRUE,
                             sort = FALSE) %>%
      arrange(factor(Original_ID, levels = anno$Original_ID)) %>%
      `rownames<-`(.$Original_ID)
  }

  ne <- .ndspaExperiment(val.All, colData = anno, rowData = probes)
  S4Vectors::metadata(ne) <- exp_type

  return(ne)
}

######## Slicing out annotation #########################################
.cut_anno <- function(df) {
  anno <- df %>%
    dplyr::slice(seq_len(which(df$`Segment displayed name` == "#Probe Group" |
                df$`Segment displayed name` == "#Target Group")) - 1) %>%
    tibble::column_to_rownames(var = "Segment displayed name") %>%
    dplyr::select(4:dim(.)[2]) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    t(.) %>%
    as.data.frame(drop = FALSE) %>%
    dplyr::mutate(Original_ID = rownames(.)) %>%
    dplyr::mutate(ID = gsub(" ", "_", paste(.$`Scan name`, .$`ROI (label)`,
                                            .$`Segment tags`, sep = "_")))

  return(anno)
}

## Slicing out Probes
.cut_probes <- function(df) {
  probes <- df %>%
    dplyr::slice(-(seq_len(which(df$`Segment displayed name` == "#Probe Group" |
    df$`Segment displayed name` == "#Target Group")) - 1)) %>%
    dplyr::select((seq_len(4))) %>%
    `colnames<-`(.[1, ]) %>%
    dplyr::slice(-1)

  return(probes)
}

## Slicing out the matrix from Data
.cut_vals <- function(df, anno, probes) {
  val <- df %>%
    dplyr::slice(-(seq_len(which(df$`Segment displayed name` == "#Probe Group" |
  df$`Segment displayed name` == "#Target Group")))) %>%
    select(-(seq_len(4))) %>%
    `colnames<-`(anno$ID) %>%
    apply(2, as.numeric) %>%
    data.matrix() %>%
    `rownames<-`(probes$`ProbeName (display name)`)

  return(val)
}

.readfn <- function(type, file) {
  if (type == "tsv") {
    suppressMessages(readr::read_tsv(file = file, col_names = TRUE))
  } else if (type == "xlsx") {
    df <- suppressMessages(readxl::read_excel(path = file, col_names = TRUE))
  } else if (type == "csv") {
    df <- suppressMessages(readr::read_csv(file = file))
  } else {
    stop("The selected file type is not available.
         Input must be in csv,tsv or xlsx format.")
  }
}

if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Original_ID"))


# styler:::style_active_file()
