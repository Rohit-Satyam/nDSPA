#' nDSPA
#'
#' Enables to perform group-wise comparison and estimation of statistically significant genes.
#'
#' The ndspaStats function is let user identify Differentially expresed genes for a given segment and groups.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @importFrom multcomp glht adjusted
#' @importFrom lmerTest lmer contest
#' @importFrom stats p.adjust
#' @importFrom lme4 isSingular
#' @export
#'
#' @param seobj S4Vector Provide ndspaExperiment object.
#' @param group.df data.frame Provide a dataframe object containing group information. Default NULL.In case user manually adds the group information, keep this option NULL.
#' @param group.by character Column name common to ndspaExperiment colData and group dataframe. This column will be used to add group information to the colData.
#' @param order vector A character vector to define the order in which group1 and group2 should be considered.
#' @param segment.tag character Segment tags to be used for performing statistical analysis.
#' @param use.assay character Provide the assay name to be used for statistical analysis. Default is "counts" assay (Raw Counts).
#' @param lowexprfilt logical To filter out low expressing genes based on min.expr (minimum Signal-to-noise ratio expression; Threshold=1) and min.sm (number of samples required to meet the threshold; Default: 3).
#' Default: "FALSE"
#' @param test character Select which statistical test should be used. Present you can use "z.test" or "t.test".
#' @param FDR logical Perform p.value correction. Default "TRUE". (This argument is not yet implemented)
#' @param p.value numeric Adjusted p-value cutoff to be used to identify significant Differentially Expressed genes. Default: "0.05"
#' @param fc numeric Specify fold change value required for a gene to qualify as Differentially Expressed. Default: "1.5"
#'
#' @return stats.df Returns a data.frame object with lFC and p.adjusted values for significantly Differentially Expressed Genes
#' @examples
#' \dontrun{
#' # To plot all QC graphs
#' dsp.group = read.delim("02-3.dsp_group.sim.tsv",
#' stringsAsFactors = F, header = T)
#' results <- ndspaStats(test,group.df=dsp.group,group.by="Scan_ID",
#' order=c("NR","R"),segment.tag = "CD45+",test="t.test")
#'
#'
#' }


ndspaStats <- function(seobj, group.df = NULL, group.by = NULL, order = NULL, segment.tag = "CD45+", use.assay = "counts", lowexprfilt = FALSE, test = c("z.test", "t.test"), FDR = TRUE, p.value = 0.05, fc = 1.5) {
  anno <- data.frame(colData(seobj), check.names = FALSE)
  temp <- ifelse(is.null(group.df) & is.null(group.by), anno <- anno, anno <- dplyr::left_join(anno, group.df, by = group.by))
  if ("group" %in% tolower(colnames(anno))) {
    cat("Groups Identified\n")
  } else {
    stop("Provide Sample Group Information\n")
  }
  val_all <- assay(seobj, use.assay)
  probes <- data.frame(rowData(seobj), check.names = FALSE)
  endo <- val_all[which(probes$`#CodeClass` == "Endogenous"), ]
  ercc <- val_all[which(probes$`#CodeClass` == "Positive" & probes$`#Analyte type` == "SpikeIn"), , drop = F]
  negprob <- val_all[which(probes$`#CodeClass` == "Negative" & (probes$`#Analyte type` == "RNA" | probes$`#Analyte type` == "Protein")), ]
  negprob <- data.frame(NegPrb.geomean = apply(negprob, 2, .gmean))
  for (i in 1:ncol(endo)) {
    endo[, i] <- (endo[, i] * anno$scalefactor[i]) / negprob$NegPrb.geomean[i]
  }
  endo <- data.frame(t(endo), check.names = F)
  endo$Original_ID <- rownames(endo)
  dsp.df <- merge(anno, endo, by = "Original_ID")
  endo <- dplyr::select(endo, -"Original_ID")
  genes <- colnames(endo)
  ## between patient
  segments <- sort(unique(dsp.df$`Segment tags`))
  if (is.null(order)) {
    group1 <- unique(dsp.df$group)[1]
    group2 <- unique(dsp.df$group)[2]
    my.seg <- segment.tag
    cat(paste0("Computing Statistics for ", my.seg, ".", group1, "_vs_", group2, "\n"))
  } else {
    group1 <- order[1]
    group2 <- order[2]
    my.seg <- segment.tag
    cat(paste0("Computing Statistics for ", my.seg, ".", group1, "_vs_", group2, "\n"))
  }
  data.plot <- dsp.df[dsp.df$`Segment tags` %in% my.seg, , drop = F]
  cat(paste0("The minimum gene expresion value was found to be: ", "\n", min(dsp.df[, genes]), "\n"))
  cat(paste0("The maximum gene expresion value was found to be: ", "\n", max(dsp.df[, genes]), "\n"))
  if (isTRUE(lowexprfilt)) {
    min.expr <- 1
    min.sm <- 3
  }

  ## Run Stats

  data.stats <- NULL
  for (i in 1:length(genes)) {
    my.gene <- genes[i]
    print(my.gene)
    my.df <- NULL
    my.stats <- NULL
    my.df <- data.plot[, c("Scan_ID", "Segment tags", "ROI_ID", "group", my.gene)]
    colnames(my.df)[5] <- "Gene"
    if (isTRUE(lowexprfilt)) {
      if (sum(my.df$Gene > min.expr, na.rm = T) >= min.sm) {
        ## keep this gene!
      } else {
        print(paste0(my.gene, " was filtered out!"))
        next
      }
    }
    ## ----------------------------------------------------------------

    ## set levels for factor in the correct order ... group1 vs group2
    my.df$group <- factor(my.df$group, levels = c(group1, group2))

    ## ----------------------------------------------------------------

    ## log transform if needed
    my.df$Gene <- log2(my.df$Gene)

    ## ----------------------------------------------------------------
    ## calculate fold change
    mean1 <- NULL
    mean2 <- NULL
    log2fc <- NULL
    ## take the avg within each sample across ROIs, then take avg across all samples
    for (my.sm in sort(unique(my.df$Scan_ID))) {
      mean1 <- c(
        mean1,
        mean(my.df$Gene[my.df$group == group1 & my.df$Scan_ID == my.sm], na.rm = T)
      )
      mean2 <- c(
        mean2,
        mean(my.df$Gene[my.df$group == group2 & my.df$Scan_ID == my.sm], na.rm = T)
      )
    }

    mean1 <- mean(mean1, na.rm = T)
    mean2 <- mean(mean2, na.rm = T)
    ## since data were log transformed, log2fc is just mean1 - mean2 (and not ratio!!)
    log2fc <- mean1 - mean2
    ## ----------------------------------------------------------------

    # fit a linear mixed effects model
    ## random effect: Scan_ID
    ## fixed effect: group

    lmer.out <- lmerTest::lmer(Gene ~ 0 + group + (1 | Scan_ID), data = my.df, REML = T)
    # parameter estimates
    summary(lmer.out)
    # set up contrasts: group1 vs group2
    contrast.matrix <- rbind(
      group1.vs.group2 = c(-1, 1)
    )

    # run contrast

    if (test == "z.test") {
      # z-test if sample size is large!
      comp.z <- multcomp::glht(lmer.out, contrast.matrix)
      # summary(comp.z, test = adjusted("none"))
      summary(comp.z, test = adjusted("BH"))
      comp.stats <- comp.z
    }
    if (test == "t.test") {
      # t-test if sample size is small!

      comp.t <- contest(lmer.out,
        L = contrast.matrix, joint = FALSE,
        ddf = "Satterthwaite",
        check_estimability = TRUE
      )
      comp.t$p.adj <- p.adjust(comp.t$`Pr(>|t|)`, method = "BH")
      comp.t

      comp.stats <- comp.t
    }
    ## ----------------------------------------------------------------

    data.stats <- rbind(
      data.stats,
      data.frame(
        Gene = my.gene,
        segment = my.seg,
        Comp = paste0(group1, "_vs_", group2),
        Group1.smtotal = sum(my.df$group == group1),
        Group2.smtotal = sum(my.df$group == group2),
        Mean1 = mean1,
        Mean2 = mean2,
        logFC = log2fc,
        comp.stats,
        is.singular = isSingular(lmer.out),
        stringsAsFactors = F
      )
    )

    rm(lmer.out)
    rm(comp.t)
    rm(comp.z)
    rm(comp.stats)
  }

  names(data.stats)[names(data.stats) == "Pr...t.."] <- "p.value"
  data.stats$p.adj <- p.adjust(data.stats$p.value, method = "fdr")
  row.names(data.stats) <- 1:nrow(data.stats)
  stats.df <- data.stats[data.stats$p.adj < p.value & abs(data.stats$logFC) >= log2(fc), ]
  rownames(stats.df) <- seq_len(nrow(stats.df))
  return(stats.df)
}


if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Scan_ID", "#CodeClass", "#Analyte type", "ProbeName (display name)", "Segment tags", "ROI_ID", "group", "len"))
#styler:::style_active_file()
