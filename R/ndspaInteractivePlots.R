#' nDSPA
#'
#' Enables visualisation of QC Plots
#'
#' The ndspaInteractivePlots function is let user generate QC plots to guage the data quality.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import dplyr
#' @import factoextra
#' @import plotly
#' @import ggplot2
#' @import shinythemes
#' @import shiny
#' @importFrom stats prcomp
#' @importFrom GGally ggpairs
#' @export
#'
#' @param seobj S4Vector Provide ndspaExperiment object.
#' @param use.assay character Provide the assay name to be used for QC plots. Default "counts".
#'
#' @return A Shiny app object is returned for interactive data exploration of ndspaExperiment object.
#' @examples
#' \dontrun{
#' # To plot all QC graphs
#' ndspaInteractivePlots(test)
#'
#' # To plot all QC graphs using erccScaled assay
#' ndspaInteractivePlots(test,use.assay ="erccScaled")
#'
#' }


ndspaInteractivePlots <- function(seobj, use.assay = "counts") {
  anno <- data.frame(SummarizedExperiment::colData(seobj), check.names = FALSE)
  val_all <- assay(seobj, assay = use.assay)
  probe <- data.frame(rowData(seobj), check.names = FALSE)
  Endo_probes <- probe$`ProbeName (display name)`[probe$`#CodeClass` == "Endogenous"]
  val_Endo <- val_all[rownames(val_all) %in% Endo_probes, ]

  shinyApp(
    ui = fluidPage(
      theme = shinytheme("flatly"),
      # Title of the app
      titlePanel("Data Plots"),
      # sidebarPanel("This Function makes 7 Data Plots",width=2),
      mainPanel(
        tabsetPanel(
          tabPanel(h5("PCA Probes", style = "font-family:'Trebuchet MS'"), plotlyOutput("pca_var", width = "100%", height = "600px")),
          tabPanel(h5("PCA Samples", style = "font-family:'Trebuchet MS'"), plotlyOutput("pca_ind", width = "100%", height = "600px")),
          tabPanel(
            h5("Density", style = "font-family:'Trebuchet MS'"), tags$br(),
            sidebarPanel(selectInput("category", "Select category to group plots", choices = as.list(c("Scan_ID","Segment tags")), selected = "Scan_ID")),
            mainPanel(plotlyOutput("density", width = "100%", height = "600px"))
          ),
          tabPanel(
            h5("Heatmap", style = "font-family:'Trebuchet MS'"), tags$br(),
            sidebarPanel(
              sliderInput("height1", h3("Height"),
                min = 0, max = 5000, value = 1000
              ),
              sliderInput("width1", h3("Width"),
                min = 0, max = 5000, value = 1500
              ),
              sliderInput("fontsize", h3("Font Size"),
                min = 0, max = 100, value = 15
              ),
              textInput("color", h3("Color Scheme"),
                value = "RdBu"
              )
            ),
            mainPanel(plotlyOutput("mainHeat"), verbatimTextOutput("code1"))
          ),
          tabPanel(h5("BG vs HK", style = "font-family:'Trebuchet MS'"), plotlyOutput("IsovHK_plt"), width = "100%", height = "800px"),
          tabPanel(h5("HK Corr", style = "font-family:'Trebuchet MS'"), plotlyOutput("HK_Sct_mt"), width = "100%", height = "800px"),
          tabPanel(h5("SNR Levels", style = "font-family:'Trebuchet MS'"), plotlyOutput("SNR_boxplot", width = "100%", height = "700px"))
        )
      )
    ),
    server = function(input, output, session) {
      ################# TAB 2 ############################################
      output$pca_ind <- renderPlotly({
        val_Endo <- val_Endo
        anno <- anno
        p <- data.matrix(val_Endo) %>%
          t() %>% apply(2, log) %>%
          stats::prcomp() %>%
          factoextra::fviz_pca_ind(label = "ind", title = "Samples", geom = c("point"), habillage = anno$Scan_ID, addEllipses = TRUE, ellipse.level = 0.95) + ggplot2::theme_minimal()
        pp <- plotly::ggplotly(p) %>% plotly::plotly_build()
      })
      ################# TAB 2 ############################################
      output$pca_var <- renderPlotly({
        val_Endo <- val_Endo
        p <- data.matrix(val_Endo) %>%
          t() %>% apply(2, log) %>%
          prcomp() %>%
          fviz_pca_var(label = "var", title = "Probes", geom = c("point", "text"), col.var = "contrib") + theme_minimal()
        plotly::ggplotly(p) %>% plotly::plotly_build()
      })
      ################# TAB 3 ###############################################
      output$density <- renderPlotly({
        anno$ID <- anno$Original_ID
        p <- data.matrix(val_Endo) %>%
          log2() %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "row.name") %>%
          tidyr::gather(key = "ID", value = "expr", -row.name) %>%
          left_join(anno, by = "ID") %>%
          ggplot() +
          geom_density(aes(x = expr, color = ID), show.legend = FALSE) +
          ggplot2::theme_minimal() +
          ggplot2::facet_wrap(paste0("`",input$category,"`"))
        plotly::ggplotly(p) %>% plotly::plotly_build()
      })
      ################# TAB 4 ######################################################
      output$mainHeat <- renderPlotly({
        val_Endo %>%
          data.matrix() %>%
          log2() %>%
          t() %>%
          heatmaply::heatmaply(
            main = "Heatmap", colorbar_xanchor = "left", colorbar_yanchor = "middle", colorbar_ypos = 0.5, row_dend_left = F,
            row_side_colors = anno[, "Scan_ID", drop = FALSE],
            column_text_angle = 90,
            fontsize_row = input$fontsize,
            fontsize_col = input$fontsize,
            plot_method = "plotly",
            scale = "column",
            colors = input$color, side_color_colorbar_len = 0.2
          ) %>%
          layout(width = input$width1, height = input$height1)
      })

      ################### TAB 5 ######################################################
      output$IsovHK_plt <- renderPlotly({
        HK <- probe %>%
          filter(`#Analyte type` == "RNA" & `#CodeClass` == "Control") %>%
          select(`ProbeName (display name)`) %>%
          unlist() %>%
          unname()
        Iso <- probe %>%
          filter(`#Analyte type` == "RNA" & `#CodeClass` == "Negative") %>%
          select(`ProbeName (display name)`) %>%
          unlist() %>%
          unname()

        Isos <- val_all[Iso, ] %>%
          apply(2, as.numeric) %>%
          apply(2, log) %>%
          apply(2, mean) %>%
          exp()

        HKs <- val_all[HK, ] %>%
          apply(2, as.numeric) %>%
          apply(2, log) %>%
          apply(2, mean) %>%
          exp()

        df <- data.frame(Isos, HKs)
        p <- ggplot(df) +
          geom_point(aes(x = HKs, y = Isos)) +
          xlab("HouseKeeping (Geomean)") +
          ylab("Background Negative (Geomean)")

        ggplotly(p) %>% layout(autosize = T)
      })
      ###################### TAB 6 #############################################
      output$HK_Sct_mt <- renderPlotly({
        # Conditional switch for val_all/val_all_QC needed
        HK <- probe %>%
          filter(`#Analyte type` == "RNA" & `#CodeClass` == "Control") %>%
          select(`ProbeName (display name)`) %>%
          unlist() %>%
          unname()


        HKs <- val_all[HK, ] %>%
          t() %>%
          as.data.frame()
        p <- GGally::ggpairs(HKs)
        ggplotly(p)
      })
      ############################ TAB 7 ############################
      output$SNR_boxplot <- renderPlotly({
        anno$ID <- anno$Original_ID
        d <- dspSNR(seobj, use.assay = use.assay)
        val_all_SNR <- assay(d, use.assay)
        p <- val_all_SNR[Endo_probes, ] %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "probe") %>%
          tidyr::gather(key = "ID", value = "expr", -probe) %>%
          mutate(expr = as.numeric(as.character(expr))) %>%
          left_join(anno, by = "ID") %>%
          ggplot(aes(x = ROI_ID, y = log2(expr), fill = Scan_ID)) +
          geom_boxplot() +
          theme_bw() +
          labs(subtitle = "SNR, log2(SNR)", title = "Comparison of ROIs per Scan", x = "ROI") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x = element_blank()) +
          facet_grid(Scan_ID ~ .)
        ggplotly(p) %>% plotly_build()
      })
    }
  )
}
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "Scan_ID", "#CodeClass", "#Analyte type", "ID", "ProbeName (display name)", "Segment tags", "ROI_ID", "group", "len", "row.name"))
#styler:::style_active_file()
