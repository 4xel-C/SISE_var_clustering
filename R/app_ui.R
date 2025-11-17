#' Shiny UI for HClustVar
#'
#' @return A Shiny UI object
#' @noRd
app_ui <- function() {
  shiny::fluidPage(

    # Titre de l'application
    shiny::titlePanel(
      shiny::div(
        style = "text-align: center; padding: 20px;",
        shiny::h1("Variable Clustering",
                  style = "color: #2C3E50; font-weight: bold;"),
        shiny::h4("Interactive application for variable clustering",
                  style = "color: #7F8C8D;")
      )
    ),

    # CSS personnalisé
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .btn-primary {
          background-color: #3498DB;
          border-color: #2980B9;
        }
        .btn-primary:hover {
          background-color: #2980B9;
        }
        .well {
          background-color: #ECF0F1;
          border: 1px solid #BDC3C7;
        }
        .progress-bar {
          background-color: #2ECC71;
        }
      "))
    ),

    # ============================================================
    # SIDEBAR - Parameters
    # ============================================================
    shiny::sidebarLayout(

      shiny::sidebarPanel(
        width = 3,

        # ------------------------------------------------------------
        # Data imports
        # ------------------------------------------------------------
        shiny::wellPanel(
          shiny::h4("📁 Data Import", style = "color: #2C3E50;"),

          shiny::fileInput(
            "file_input",
            "Upload Excel/CSV file",
            accept = c(".xlsx", ".xls", ".csv"),
            placeholder = "No file selected"
          ),

          shiny::conditionalPanel(
            condition = "output.file_uploaded",
            shiny::checkboxInput(
              "header",
              "First row as column names",
              value = TRUE
            ),

            shiny::numericInput(
              "sheet_index",
              "Excel sheet number",
              value = 1,
              min = 1,
              step = 1
            ),

            shiny::br(),

            # Button to open variable configuration modal
            shiny::actionButton(
              "open_var_config",
              "⚙️ Configure Variables",
              class = "btn-info btn-block",
              style = "font-weight: bold;"
            )
          )
        ),

        # ------------------------------------------------------------
        # Clustering algorithm selection
        # ------------------------------------------------------------
        shiny::wellPanel(
          shiny::h4("⚙️ Clustering Parameters", style = "color: #2C3E50;"),

          shiny::selectInput(
            "algorithm",
            "Algorithm",
            choices = c(
              "CAH Variables clustering" = "hclust",
              "K-means  Variables clustering" = "kmeans"
            ),
            selected = "hclust"
          ),


          # ------------------------------------------------------------
          # Type selection
          # ------------------------------------------------------------
          shiny::conditionalPanel(
            condition = "input.algorithm == 'hclust'",
            shiny::selectInput(
              "vartype",
              "Variable type",
              choices = c(
                "Auto-detect" = "auto",
                "Quantitative" = "quant",
                "Qualitative" = "qual",
                "Mixed" = "mixed"
              ),
              selected = "auto"
            )
          ),


          shiny::conditionalPanel(
            condition = "input.algorithm == 'kmeans'",
            shiny::selectInput(
              "vartype_kmeans",
              "Variable type",
              choices = c(
                "Quantitative" = "quant"
              ),
              selected = "quant"
            )
          ),

          # ------------------------------------------------------------
          # Descriptive text for selection
          # ------------------------------------------------------------
          shiny::conditionalPanel(
            condition = "input.vartype == 'auto' && input.algorithm == 'hclust'",
            shiny::div(
              style = "background-color: #E8F4F8; padding: 10px; border-radius: 5px; margin-top: 10px; border-left: 3px solid #3498DB;",
              shiny::tags$small(
                shiny::icon("info-circle"),
                " The algorithm will automatically detect the type of variables in your dataset."
              )
            )
          ),

          # Descriptive text for quantitative
          shiny::conditionalPanel(
            condition = "(input.vartype == 'quant' && input.algorithm == 'hclust') || (input.vartype_kmeans == 'quant' && input.algorithm == 'kmeans') ",
            shiny::div(
              style = "background-color: #E8F5E9; padding: 10px; border-radius: 5px; margin-top: 10px; border-left: 3px solid #2ECC71;",
              shiny::tags$small(
                shiny::icon("chart-line"),
                shiny::strong(" Quantitative variables:"),
                "Compute clustering on quantitatives data only"
              )
            )
          ),

          # Descriptive text for qualitative
          shiny::conditionalPanel(
            condition = "input.vartype == 'qual' && input.algorithm == 'hclust'",
            shiny::div(
              style = "background-color: #FFF3E0; padding: 10px; border-radius: 5px; margin-top: 10px; border-left: 3px solid #F39C12;",
              shiny::tags$small(
                shiny::icon("tags"),
                shiny::strong(" Qualitative variables:"),
                "Compute clustering on categorical data only"
              )
            )
          ),

          # Descriptive text for mixed
          shiny::conditionalPanel(
            condition = "input.vartype == 'mixed' && input.algorithm == 'hclust'",
            shiny::div(
              style = "background-color: #F3E5F5; padding: 10px; border-radius: 5px; margin-top: 10px; border-left: 3px solid #9B59B6;",
              shiny::tags$small(
                shiny::icon("layer-group"),
                shiny::strong(" Mixed variables:"),
                "Will select both quantitative and qualitative variables. If mixed type, proceed a discretization of quantitatives."
              )
            )
          ),

          # ------------------------------------------------------------
          # Parameters for Correlation
          # ------------------------------------------------------------
          shiny::conditionalPanel(
            condition = "(input.algorithm == 'hclust' && input.vartype == 'quant') || input.algorithm == 'kmeans'",
            shiny::selectInput(
              "dist_metric",
              "Distance metric",
              choices = c(
                "Squared correlation (R²)" = "rsquare",
                "Correlation (r)" = "r"
              ),
              selected = "rsquare"
            )
          ),


          # ------------------------------------------------------------
          # Parameters for CAH
          # ------------------------------------------------------------
          shiny::conditionalPanel(
            condition = "input.algorithm == 'hclust'",
            shiny::selectInput(
              "cah_method",
              "Clustering method",
              choices = c(
                "Ward D" = "ward.D",
                "Ward D2" = "ward.D2",
                "Complete" = "complete",
                "Average" = "average",
                "Single" = "single",
                "Mcquitty" = "mcquitty",
                "Median" = "median",
                "Centroid" = "centroid"
              ),
              selected = "ward.D"
            )
          ),


          # ------------------------------------------------------------
          # Parameters for Kmeans
          # ------------------------------------------------------------
          shiny::conditionalPanel(
            condition = "input.algorithm == 'kmeans'",
            shiny::numericInput(
              "km_max_iter",
              "K-means max iterations",
              value = 100,
              min = 1,
              step = 1
            ),
            shiny::numericInput(
              "km_n_init",
              "K-means n_init (restarts)",
              value = 10,
              min = 1,
              step = 1
            ),
            shiny::numericInput(
              "km_random_state",
              "Random seed (0 = random)",
              value = 0,
              min = 0,
              step = 1
            )
          ),

          # ------------------------------------------------------------
          # Parameters for Number of clusters
          # ------------------------------------------------------------

          # --- Common parameter: number of clusters ---
          shiny::numericInput(
            "n_clusters",
            "Number of clusters",
            value = 3,
            min = 2,
            max = 10,
            step = 1
          ),


          shiny::actionButton(
            "run_clustering",
            "🚀 Run Clustering",
            class = "btn-primary btn-block",
            style = "margin-top: 15px; font-weight: bold;"
          )
        ),

        # ------------------------------------------------------------
        # Data export
        # ------------------------------------------------------------
        shiny::wellPanel(
          shiny::h4("💾 Export Results", style = "color: #2C3E50;"),

          shiny::conditionalPanel(
            condition = "output.clustering_done",
            shiny::downloadButton(
              "download_results",
              "Download Summary (CSV)",
              class = "btn-success btn-block"
            )
          )
        )
      ),

      # ============================================================
      # MAIN PANEL - Dynamic results
      # ============================================================
      shiny::mainPanel(
        width = 9,

        # Les tabs seront générés dynamiquement côté serveur
        shiny::uiOutput("dynamic_tabs")
      )
    )
  )
}
