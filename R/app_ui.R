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
        shiny::h1("Hierarchical Variable Clustering",
                  style = "color: #2C3E50; font-weight: bold;"),
        shiny::h4("Interactive Analysis with HClustVar",
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
    # SIDEBAR - Paramètres
    # ============================================================
    shiny::sidebarLayout(

      shiny::sidebarPanel(
        width = 3,

        # --- SECTION 1: Import de données ---
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
            )
          )
        ),

        # --- SECTION 2: Paramètres de clustering ---
        shiny::wellPanel(
          shiny::h4("⚙️ Clustering Parameters", style = "color: #2C3E50;"),

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
          ),

          shiny::conditionalPanel(
            condition = "input.vartype == 'quant'",
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
          ),

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

        # --- SECTION 3: Export ---
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
      # MAIN PANEL - Résultats
      # ============================================================
      shiny::mainPanel(
        width = 9,

        shiny::tabsetPanel(
          id = "main_tabs",
          type = "tabs",

          # --- TAB 1: Data Preview ---
          shiny::tabPanel(
            "📊 Data Preview",
            shiny::br(),
            shiny::uiOutput("data_info"),
            shiny::br(),
            DT::DTOutput("data_table")
          ),

          # --- TAB 2: Dendrogram ---
          shiny::tabPanel(
            "🌳 Dendrogram",
            shiny::br(),
            shiny::plotOutput("dendrogram_plot", height = "600px")
          ),

          # --- TAB 3: Aggregation Levels ---
          shiny::tabPanel(
            "📈 Elbow Method",
            shiny::br(),
            shiny::plotOutput("agg_levels_plot", height = "600px")
          ),

          # --- TAB 4: Results Summary ---
          shiny::tabPanel(
            "📋 Summary",
            shiny::br(),
            shiny::h3("Cluster Summary"),
            DT::DTOutput("cluster_summary_table"),
            shiny::br(),
            shiny::h3("Variable Assignments"),
            DT::DTOutput("cluster_members_table")
          ),

          # --- TAB 5: Silhouette ---
          shiny::tabPanel(
            "📊 Silhouette Plot",
            shiny::br(),
            shiny::plotOutput("silhouette_plot", height = "700px")
          ),

          # --- TAB 6: MDS Projection ---
          shiny::tabPanel(
            "🗺️ MDS Projection",
            shiny::br(),
            shiny::plotOutput("mds_plot", height = "600px")
          )
        )
      )
    )
  )
}
