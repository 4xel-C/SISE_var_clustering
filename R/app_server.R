#' Shiny Server for HClustVar
#'
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#'
#' @return Server function
#' @noRd
app_server <- function(input, output, session) {

  # ===========================================================
  # REACTIVE VALUES
  # ===========================================================

  rv <- shiny::reactiveValues(
    data = NULL,
    model = NULL,
    model_type = NULL,
    clustering_done = FALSE,
    prediction = NULL,
    prediction_proximities = NULL
  )

  # Reactive value to store variable configurations
  variable_config <- shiny::reactiveVal(NULL)


  # ===========================================================
  # DYNAMIC TAB GENERATION
  # ===========================================================

  output$dynamic_tabs <- shiny::renderUI({

    algorithm <- input$algorithm
    if (is.null(algorithm)) algorithm <- "hclust"

    # ------------------------------------------------------------
    # Common tabs
    # ------------------------------------------------------------
    common_tabs <- list(

      # -------------
      # Data preview
      # -------------
      shiny::tabPanel(
        "📊 Data Preview",
        shiny::br(),
        shiny::uiOutput("data_info"),
        shiny::br(),
        DT::DTOutput("data_table")
      ),

      # -------------
      # Variable configuration
      # -------------
      shiny::tabPanel(
        "🔍 Variable Configuration",
        shiny::br(),

        # Instructions and button
        shiny::wellPanel(
          style = "background-color: #E3F2FD; border-left: 4px solid #2196F3;",
          shiny::fluidRow(
            shiny::column(
              8,
              shiny::h4("📋 Variable Configuration", style = "color: #1976D2; margin-top: 0;"),
              shiny::tags$p("Click the button to configure which variables to use in clustering and their types.")
            ),
            shiny::column(
              4,
              shiny::actionButton(
                "reopen_var_config",
                "⚙️ Configure Variables",
                class = "btn-primary btn-lg btn-block",
                style = "margin-top: 10px; font-weight: bold;"
              )
            )
          )
        ),

        shiny::br(),

        # Summary of current configuration
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::wellPanel(
              style = "background-color: #E8F5E9; min-height: 300px;",
              shiny::h4("✅ Active Variables", style = "color: #2ECC71;"),
              shiny::tags$p("Variables used in clustering analysis"),
              shiny::hr(),
              shiny::uiOutput("active_vars_display")
            )
          ),
          shiny::column(
            4,
            shiny::wellPanel(
              style = "background-color: #FFF3E0; min-height: 300px;",
              shiny::h4("ℹ️ Illustrative Variables", style = "color: #F39C12;"),
              shiny::tags$p("Variables for interpretation only"),
              shiny::hr(),
              shiny::uiOutput("illustrative_vars_display")
            )
          ),
          shiny::column(
            4,
            shiny::wellPanel(
              style = "background-color: #FFEBEE; min-height: 300px;",
              shiny::h4("❌ Excluded Variables", style = "color: #E74C3C;"),
              shiny::tags$p("Variables completely ignored"),
              shiny::hr(),
              shiny::uiOutput("excluded_vars_display")
            )
          )
        ),

        shiny::br(),

        # Detailed configuration table (read-only)
        shiny::wellPanel(
          shiny::h4("📊 Complete Configuration"),
          DT::DTOutput("config_summary_table")
        )
      ),

      # -------------
      # Summary
      # -------------
      shiny::tabPanel(
        "📋 Summary",

        # --------------
        # Metrics cards
        # --------------
        shiny::fluidRow(
          # Average Silhouette Card
          shiny::column(
            6,
            shiny::wellPanel(
              style = "background-color: #E3F2FD; border-left: 5px solid #2196F3; min-height: 120px;",
              shiny::div(
                style = "text-align: center;",
                shiny::h4("📊 Average Silhouette", style = "color: #1976D2; margin-top: 0;"),
                shiny::uiOutput("avg_silhouette_display")
              )
            )
          ),

          # Total Variance Explained Card
          shiny::column(
            6,
            shiny::wellPanel(
              style = "background-color: #E8F5E9; border-left: 5px solid #4CAF50; min-height: 120px;",
              shiny::div(
                style = "text-align: center;",
                shiny::h4("📈 Total Variance Explained", style = "color: #2E7D32; margin-top: 0;"),
                shiny::uiOutput("total_var_explained_display")
              )
            )
          )
        ),

        # --------------
        # Cluster summary
        # --------------
        shiny::wellPanel(
          shiny::h3("Cluster Summary"),
          DT::DTOutput("cluster_summary_table")
        ),

        shiny::br(),

        # --------------
        # Centroid correlation
        # --------------
        shiny::wellPanel(
          shiny::h3("Centroids Correlations"),
          shiny::tags$p("Correlation matrix between cluster centroids"),
          DT::DTOutput("centroids_correlation_table")
        ),

        shiny::br(),

        # --------------
        # Variable assignements
        # --------------
        shiny::wellPanel(
          shiny::h3("Variable Assignments"),
          DT::DTOutput("cluster_members_table")
        ),

        # --------------
        # All clusters R² matrix  <-- NOUVEAU
        # --------------
        shiny::wellPanel(
          shiny::h3("R² Matrix: All Variables vs All Clusters"),
          shiny::tags$p(
            "Complete correlation matrix showing R² values between each variable and all cluster centroids.",
            shiny::tags$br(),
            "Higher values indicate stronger association with that cluster."
          ),
          DT::DTOutput("all_clusters_R2_table")
        )
      )
    )



    # ------------------------------------------------------------
    # Prediction tab
    # ------------------------------------------------------------
    predict_tab <- list(
      shiny::tabPanel(
        "➡️ Predict",
        shiny::br(),

        # Limit the width of the table.
        shiny::fluidRow(
          shiny::column(
            10,
            offset = 1,

            shiny::h3("Cluster prediction on illustrative variables"),
            shiny::br(),

            # Message if no prediction
            if (is.null(rv$prediction)) {
              shiny::wellPanel(
                style = "background-color: #FFF3E0; border-left: 4px solid #F39C12;",
                shiny::h4(
                  "ℹ️ No predictions available",
                  style = "color: #F39C12; margin-top: 0;"
                ),
                shiny::tags$p("Configure Illustrative Variables and run a clustering to make your prediction.")
              )
            },

            # Prediction labels table
            shiny::wellPanel(
              shiny::h4("📌 Predicted Clusters"),
              shiny::tags$p("Cluster assignment for each illustrative variable"),
              DT::DTOutput("prediction_result")
            ),

            shiny::br(),

            # Proximities table
            shiny::wellPanel(
              shiny::h4("📊 Proximity Matrix"),
              shiny::tags$p(
                "Correlation/proximity values between each illustrative variable and all cluster centroids.",
                shiny::tags$br(),
                "Higher values indicate stronger association with that cluster."
              ),
              DT::DTOutput("prediction_proximities_table")
            )
          )
        )
      )
    )

    # ------------------------------------------------------------
    # HClust tabs
    # ------------------------------------------------------------
    hclust_tabs <- list(

      # -------------
      # Silhouette plot
      # -------------
      shiny::tabPanel(
        "📊 Silhouette Plot",
        shiny::br(),
        shiny::plotOutput("silhouette_plot", height = "700px")
      ),

      shiny::tabPanel(
        "🌳 Dendrogram",
        shiny::br(),
        shiny::plotOutput("dendrogram_plot", height = "600px")
      ),

      shiny::tabPanel(
        "📈 Elbow Method",
        shiny::br(),
        shiny::plotOutput("agg_levels_plot", height = "600px")
      ),

      shiny::tabPanel(
        "🗺️ Factorial Projection",
        shiny::br(),
        shiny::plotOutput("mds_plot", height = "600px")
      )
    )

    # ------------------------------------------------------------
    # Kmeans tabs
    # ------------------------------------------------------------
    kmeans_tabs <- list(
      shiny::tabPanel(
        "📉 Inertia Evolution",
        shiny::br(),
        shiny::plotOutput("inertia_plot", height = "600px")
      ),

      shiny::tabPanel(
        "🗺️ Variable Projection",
        shiny::br(),
        shiny::plotOutput("kmeans_projection_plot", height = "600px")
      ),

      shiny::tabPanel(
        "🔄 Cluster Centers",
        shiny::br(),
        shiny::h4("Cluster Centers Heatmap"),
        shiny::plotOutput("centers_heatmap", height = "600px")
      )
    )

    # TODO: ModalClust tab to be instantiated

    # ------------------------------------------------------------
    # Select the vector of tabs to display
    # ------------------------------------------------------------
    if (algorithm == "hclust") {
      all_tabs <- c(common_tabs, hclust_tabs, predict_tab)
    } else if (algorithm == "kmeans") {
      all_tabs <- c(common_tabs, kmeans_tabs, predict_tab)
    } else {
      all_tabs <- c(common_tabs, predict_tab)
    }

    # Create the tabsetpanel with the list of tabs to use
    do.call(shiny::tabsetPanel, c(list(id = "main_tabs", type = "tabs"), all_tabs))
  })

  # ===========================================================
  # DATA LOADING
  # ===========================================================

  # Load the data from the selected file.
  data_loaded <- shiny::reactive({

    req(input$file_input)

    file_path <- input$file_input$datapath
    file_ext <- tools::file_ext(input$file_input$name)

    shiny::withProgress(message = 'Loading data...', value = 0.5, {

      tryCatch({

        if (file_ext %in% c("xlsx", "xls")) {
          # Excel
          data <- readxl::read_excel(
            file_path,
            sheet = input$sheet_index,
            col_names = input$header
          )
        } else if (file_ext == "csv") {
          # CSV
          data <- read.csv(
            file_path,
            header = input$header,
            stringsAsFactors = FALSE
          )
        } else {
          stop("Unsupported file format")
        }

        # Convert into dataframe
        data <- as.data.frame(data)

        # Detect the column to convert into factor automaticly
        for (col in names(data)) {
          if (is.character(data[[col]]) ||
              (is.numeric(data[[col]]) && length(unique(data[[col]])) <= 2)) {
            data[[col]] <- as.factor(data[[col]])
          }
        }

        rv$data <- data
        rv$clustering_done <- FALSE

        return(data)

      }, error = function(e) {
        shiny::showNotification(
          paste("Error loading file:", e$message),
          type = "error",
          duration = 5
        )
        return(NULL)
      })
    })
  })


  # ===========================================================
  # DATA PREVIEW
  # ===========================================================

  output$file_uploaded <- shiny::reactive({
    return(!is.null(input$file_input))
  })
  shiny::outputOptions(output, "file_uploaded", suspendWhenHidden = FALSE)


  # ===========================================================
  # VARIABLE CONFIGURATION LISTENER
  # ===========================================================

  # Initialize variable configuration when data is loaded
  shiny::observe({
    req(data_loaded())

    data <- data_loaded()
    var_names <- names(data)

    # Create initial configuration dataframe with all variables as "Active"
    config_df <- data.frame(
      Variable = var_names,
      Type = sapply(data, function(x) {
        if (is.numeric(x)) "Quantitative" else "Qualitative"
      }),
      Role = rep("Active", length(var_names)),
      stringsAsFactors = FALSE
    )

    variable_config(config_df)
  })

  # ===========================================================
  # MODAL FOR VARIABLE CONFIG
  # ===========================================================

  # Auto-open modal when data is loaded
  shiny::observe({
    req(data_loaded())
    req(variable_config())

    # Small delay to ensure data is fully loaded
    shiny::invalidateLater(500, session)
    shiny::isolate({
      if (!is.null(input$file_input)) {
        shiny::showModal(create_variable_config_modal())
      }
    })
  }) %>%
    # execute the code only the first time the dataset is loaded.
    shiny::bindEvent(input$file_input, ignoreInit = TRUE, once = TRUE)

  # Open modal manually with button
  shiny::observeEvent(input$open_var_config, {
    req(variable_config())
    shiny::showModal(create_variable_config_modal())
  })

  # Function to create the modal dialog
  create_variable_config_modal <- function() {
    req(variable_config())

    # get the configuration dataframe of the variables.
    config_df <- variable_config()

    # Create UI elements for each variable, stored in a variable.
    variable_inputs <- lapply(1:nrow(config_df), function(i) {
      var_name <- config_df$Variable[i]
      var_type <- config_df$Type[i]
      var_role <- config_df$Role[i]

      shiny::div(
        style = "border: 1px solid #ddd; padding: 10px; margin-bottom: 10px; border-radius: 5px; background-color: #f9f9f9;",
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::tags$strong(
              style = "font-size: 14px; color: #2C3E50;",
              var_name
            )
          ),
          shiny::column(
            4,
            shiny::selectInput(
              inputId = paste0("var_type_", i),
              label = NULL,
              choices = c("Quantitative", "Qualitative"),
              selected = var_type,
              width = "100%"
            )
          ),
          shiny::column(
            4,
            shiny::selectInput(
              inputId = paste0("var_role_", i),
              label = NULL,
              choices = c("Active", "Illustrative", "Excluded"),
              selected = var_role,
              width = "100%"
            )
          )
        )
      )
    })

    shiny::modalDialog(
      title = shiny::div(
        style = "color: #2C3E50;",
        shiny::icon("cogs"),
        " Variable Configuration"
      ),
      size = "l",

      # Instructions
      shiny::wellPanel(
        style = "background-color: #E3F2FD; border-left: 4px solid #2196F3;",
        shiny::h5("📋 Instructions", style = "color: #1976D2; margin-top: 0;"),
        shiny::tags$ul(
          style = "margin-bottom: 0;",
          shiny::tags$li(shiny::strong("Type:"), " Quantitative (numeric) or Qualitative (categorical)"),
          shiny::tags$li(shiny::strong("Active:"), " Variables used in clustering"),
          shiny::tags$li(shiny::strong("Illustrative:"), " Variables for interpretation only"),
          shiny::tags$li(shiny::strong("Excluded:"), " Variables to ignore completely"),
          shiny::tags$li("⚠️ ", shiny::strong("Minimum 3 active variables required"))
        )
      ),

      # Column headers
      shiny::fluidRow(
        style = "background-color: #34495e; color: white; padding: 10px; margin-bottom: 10px; border-radius: 5px; font-weight: bold;",
        shiny::column(4, "Variable Name"),
        shiny::column(4, "Type"),
        shiny::column(4, "Role")
      ),

      # Scrollable area for variables
      shiny::div(
        style = "max-height: 400px; overflow-y: auto; padding-right: 10px;",
        variable_inputs
      ),

      # Summary footer
      shiny::br(),
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::div(
            style = "background-color: #E8F5E9; padding: 10px; border-radius: 5px; text-align: center;",
            shiny::uiOutput("modal_active_count")
          )
        ),
        shiny::column(
          4,
          shiny::div(
            style = "background-color: #FFF3E0; padding: 10px; border-radius: 5px; text-align: center;",
            shiny::uiOutput("modal_illustrative_count")
          )
        ),
        shiny::column(
          4,
          shiny::div(
            style = "background-color: #FFEBEE; padding: 10px; border-radius: 5px; text-align: center;",
            shiny::uiOutput("modal_excluded_count")
          )
        )
      ),

      footer = shiny::tagList(
        shiny::actionButton("cancel_var_config", "Cancel", class = "btn-default"),
        shiny::actionButton("save_var_config", "Save Configuration", class = "btn-primary")
      )
    )
  }

  # Dynamic counters in modal
  output$modal_active_count <- shiny::renderUI({
    req(variable_config())
    config_df <- variable_config()
    n_vars <- nrow(config_df)

    # Count current selections from inputs
    active_count <- 0
    for (i in 1:n_vars) {
      role_input <- input[[paste0("var_role_", i)]]
      if (!is.null(role_input) && role_input == "Active") {
        active_count <- active_count + 1
      }
    }

    color <- if (active_count < 3) "#E74C3C" else "#2ECC71"

    shiny::tagList(
      shiny::tags$div(
        style = paste0("color: ", color, "; font-weight: bold; font-size: 16px;"),
        paste("✅", active_count, "Active")
      )
    )
  })

  output$modal_illustrative_count <- shiny::renderUI({
    req(variable_config())
    config_df <- variable_config()
    n_vars <- nrow(config_df)

    illus_count <- 0
    for (i in 1:n_vars) {
      role_input <- input[[paste0("var_role_", i)]]
      if (!is.null(role_input) && role_input == "Illustrative") {
        illus_count <- illus_count + 1
      }
    }

    shiny::tagList(
      shiny::tags$div(
        style = "color: #F39C12; font-weight: bold; font-size: 16px;",
        paste("ℹ️", illus_count, "Illustrative")
      )
    )
  })

  output$modal_excluded_count <- shiny::renderUI({
    req(variable_config())
    config_df <- variable_config()
    n_vars <- nrow(config_df)

    excl_count <- 0
    for (i in 1:n_vars) {
      role_input <- input[[paste0("var_role_", i)]]
      if (!is.null(role_input) && role_input == "Excluded") {
        excl_count <- excl_count + 1
      }
    }

    shiny::tagList(
      shiny::tags$div(
        style = "color: #95A5A6; font-weight: bold; font-size: 16px;",
        paste("❌", excl_count, "Excluded")
      )
    )
  })

  # Save configuration
  shiny::observeEvent(input$save_var_config, {
    req(variable_config())

    config_df <- variable_config()
    n_vars <- nrow(config_df)

    # Collect all selections
    new_config <- config_df
    active_count <- 0

    for (i in 1:n_vars) {
      type_input <- input[[paste0("var_type_", i)]]
      role_input <- input[[paste0("var_role_", i)]]

      if (!is.null(type_input)) {
        new_config$Type[i] <- type_input
      }
      if (!is.null(role_input)) {
        new_config$Role[i] <- role_input
        if (role_input == "Active") active_count <- active_count + 1
      }
    }

    # Validate minimum active variables
    if (active_count < 3) {
      shiny::showNotification(
        "⚠️ Please select at least 3 active variables for clustering.",
        type = "warning",
        duration = 5
      )
      return()
    }

    # Save configuration
    variable_config(new_config)

    shiny::showNotification(
      paste("✅ Configuration saved:", active_count, "active variables"),
      type = "message",
      duration = 3
    )

    # Close modal
    shiny::removeModal()
  })

  # Cancel configuration
  shiny::observeEvent(input$cancel_var_config, {
    shiny::removeModal()
  })


  # Reopen modal from Variable Configuration tab
  shiny::observeEvent(input$reopen_var_config, {
    req(variable_config())
    shiny::showModal(create_variable_config_modal())
  })


  # ===========================================================
  # FILTER DATA, TYPES AND ILLUSTRATIVE VARIABLES
  # ===========================================================

  # Create filtered data based on variable configuration
  data_filtered <- shiny::reactive({
    req(data_loaded())
    req(variable_config())

    data <- data_loaded()

    config_df <- variable_config()

    # Filter to active variables only
    active_vars <- config_df$Variable[config_df$Role == "Active" | config_df$Role == "Illustrative"]

    if (length(active_vars) == 0) {
      return(data)  # Return all if no selection (shouldn't happen)
    }

    # Apply variable type changes
    data_filtered <- data[, active_vars, drop = FALSE]

    # Convert variable types according to configuration
    for (var in active_vars) {
      var_type <- config_df$Type[config_df$Variable == var]

      if (var_type == "Qualitative" && !is.factor(data_filtered[[var]])) {
        data_filtered[[var]] <- as.factor(data_filtered[[var]])
      } else if (var_type == "Quantitative" && !is.numeric(data_filtered[[var]])) {
        # Try to convert to numeric
        data_filtered[[var]] <- suppressWarnings(as.numeric(as.character(data_filtered[[var]])))
      }
    }

    return(data_filtered)
  })

  # ----------------------------------------------------------
  # illustrative variables
  # ----------------------------------------------------------

  # Store illustrative variables for later use
  illustrative_vars <- shiny::reactive({
    req(variable_config())
    config_df <- variable_config()

    return(config_df$Variable[config_df$Role == "Illustrative"])
  })

  # ----------------------------------------------------------
  # quantitatives variables
  # ----------------------------------------------------------

  # Store quantitatives variables for later use
  quantitative_vars <- shiny::reactive({
    req(variable_config())
    config_df <- variable_config()

    return(config_df$Variable[(config_df$Type == "Quantitative") & (config_df$Role == "Active")])
  })

  # ----------------------------------------------------------
  # qualitatives variables
  # ----------------------------------------------------------

  # Store qualitatives variables for later use
  qualitative_vars <- shiny::reactive({
    req(variable_config())
    config_df <- variable_config()

    return(config_df$Variable[(config_df$Type == "Qualitative") & (config_df$Role == "Active")])
  })


  # ===========================================================
  # DISPLAY INFORMATION ABOUT DATA
  # ===========================================================


  # ----------------------------------------------------------
  # Variable configuration tab
  # ----------------------------------------------------------

  # Display active variables
  output$active_vars_display <- shiny::renderUI({
    req(variable_config())

    config_df <- variable_config()
    active_vars <- config_df[config_df$Role == "Active", ]

    if (nrow(active_vars) == 0) {
      return(shiny::tags$p(
        style = "color: gray; font-style: italic;",
        "No active variables configured"
      ))
    }

    shiny::tagList(
      shiny::tags$p(
        style = "font-weight: bold; font-size: 16px; color: #2ECC71;",
        paste(nrow(active_vars), "variable(s)")
      ),
      shiny::tags$ul(
        style = "max-height: 200px; overflow-y: auto;",
        lapply(1:nrow(active_vars), function(i) {
          shiny::tags$li(
            shiny::tags$strong(active_vars$Variable[i]),
            " - ",
            shiny::tags$span(
              style = "color: #7F8C8D; font-size: 0.9em;",
              active_vars$Type[i]
            )
          )
        })
      )
    )
  })

  # Display illustrative variables
  output$illustrative_vars_display <- shiny::renderUI({
    req(variable_config())

    config_df <- variable_config()
    illus_vars <- config_df[config_df$Role == "Illustrative", ]

    if (nrow(illus_vars) == 0) {
      return(shiny::tags$p(
        style = "color: gray; font-style: italic;",
        "No illustrative variables configured"
      ))
    }

    shiny::tagList(
      shiny::tags$p(
        style = "font-weight: bold; font-size: 16px; color: #F39C12;",
        paste(nrow(illus_vars), "variable(s)")
      ),
      shiny::tags$ul(
        style = "max-height: 200px; overflow-y: auto;",
        lapply(1:nrow(illus_vars), function(i) {
          shiny::tags$li(
            shiny::tags$strong(illus_vars$Variable[i]),
            " - ",
            shiny::tags$span(
              style = "color: #7F8C8D; font-size: 0.9em;",
              illus_vars$Type[i]
            )
          )
        })
      )
    )
  })

  # Display excluded variables
  output$excluded_vars_display <- shiny::renderUI({
    req(variable_config())

    config_df <- variable_config()
    excl_vars <- config_df[config_df$Role == "Excluded", ]

    if (nrow(excl_vars) == 0) {
      return(shiny::tags$p(
        style = "color: gray; font-style: italic;",
        "No excluded variables"
      ))
    }

    shiny::tagList(
      shiny::tags$p(
        style = "font-weight: bold; font-size: 16px; color: #E74C3C;",
        paste(nrow(excl_vars), "variable(s)")
      ),
      shiny::tags$ul(
        style = "max-height: 200px; overflow-y: auto;",
        lapply(1:nrow(excl_vars), function(i) {
          shiny::tags$li(
            shiny::tags$strong(excl_vars$Variable[i]),
            " - ",
            shiny::tags$span(
              style = "color: #7F8C8D; font-size: 0.9em;",
              excl_vars$Type[i]
            )
          )
        })
      )
    )
  })

  # Display configuration summary table
  output$config_summary_table <- DT::renderDT({
    req(variable_config())

    config_df <- variable_config()

    DT::datatable(
      config_df,
      options = list(
        pageLength = 15,
        dom = 'ftp',
        ordering = TRUE
      ),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        'Role',
        backgroundColor = DT::styleEqual(
          c('Active', 'Illustrative', 'Excluded'),
          c('#D5F4E6', '#FFF3E0', '#FFEBEE')
        ),
        fontWeight = 'bold'
      ) %>%
      DT::formatStyle(
        'Type',
        backgroundColor = DT::styleEqual(
          c('Quantitative', 'Qualitative'),
          c('#E3F2FD', '#F3D9FF')
        )
      )
  })


  # ----------------------------------------------------------
  # Data preview tab
  # ----------------------------------------------------------

  output$data_info <- shiny::renderUI({
    req(data_loaded())

    data_full <- data_loaded()
    data <- data_filtered()

    n_rows <- nrow(data)
    n_cols <- ncol(data)
    n_numeric <- sum(sapply(data, is.numeric))
    n_factor <- sum(sapply(data, is.factor))

    # Check if variables are filtered
    is_filtered <- ncol(data_full) != ncol(data)

    shiny::wellPanel(
      style = "background-color: #D5F4E6; border-left: 5px solid #27AE60;",
      shiny::h4("✅ Data loaded successfully!", style = "color: #27AE60; margin-top: 0;"),
      shiny::tags$ul(
        shiny::tags$li(paste("Rows:", n_rows)),
        shiny::tags$li(paste("Columns:", n_cols, if(is_filtered) paste0(" (", ncol(data_full), " total)") else "")),
        shiny::tags$li(paste("Numeric variables:", n_numeric)),
        shiny::tags$li(paste("Categorical variables:", n_factor))
      ),
      if(is_filtered) {
        shiny::div(
          style = "background-color: #E3F2FD; padding: 8px; border-radius: 5px; margin-top: 10px;",
          shiny::tags$small(
            shiny::icon("filter"),
            shiny::strong(" Variable filtering active: "),
            paste(n_cols, "of", ncol(data_full), "variables selected")
          )
        )
      }
    )
  })

  # Data preview
  output$data_table <- DT::renderDT({
    req(data_filtered())

    DT::datatable(
      data_filtered(),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip'
      ),
      rownames = FALSE,
      class = 'cell-border stripe'
    )
  })


  # ===========================================================
  # CLUSTERING EXECUTION
  # ===========================================================

  shiny::observeEvent(input$run_clustering, {

    # Validate minimum number of variables
    if (is.null(rv$data)) {
      shiny::showNotification(
        "⚠️ No data loaded!",
        type = "warning",
        duration = 5
      )
      return()
    }


    req(data_filtered())
    req(variable_config())

    # Validate minimum number of variables
    if (ncol(data_filtered()) < 3) {
      shiny::showNotification(
        "⚠️ Please select at least 3 variables for clustering.",
        type = "warning",
        duration = 5
      )
      return()
    }

    shiny::withProgress(message = 'Running clustering...', value = 0, {

      # -----------------------------------------------------------
      # Error handling and data preparation
      # -----------------------------------------------------------

      tryCatch({

        # Step 1: Initialization
        shiny::incProgress(0.2, detail = "Initializing model...")

        # Get the configuration dataframe.
        config_df <- variable_config()

        # Get the filtered data (contains active + illustrative variables)
        data <- req(data_filtered())  # exit if NULL

        # Get quantitative, qualitative, illustrative and active variables
        quant_vars <- setdiff(quantitative_vars(), illustrative_vars())
        qual_vars  <- setdiff(qualitative_vars(), illustrative_vars())
        illust_vars <- illustrative_vars()
        active_vars <- c(quant_vars, qual_vars)


        # Check the selected algorithm
        if (is.null(input$algorithm) || input$algorithm == 'hclust') {

          # Choose the distance type
          dist_metric <- if (input$vartype == "quant") input$dist_metric else NULL

          # Check if there are variables for the chosen type
          if (length(quant_vars) == 0 && input$vartype == "quant") {
            shiny::showNotification("⚠️ No quantitative variables selected!", type = "warning")
            return(NULL)
          } else if (length(qual_vars) == 0 && input$vartype == "qual") {
            shiny::showNotification("⚠️ No qualitative variables selected!", type = "warning")
            return(NULL)
          }

          # Select only the variables of the chosen type
          if (input$vartype == "quant") {
            data <- data[, quant_vars, drop = FALSE]
          } else if (input$vartype == "qual") {
            data <- data[, qual_vars, drop = FALSE]
          } else {
            data <- data[, active_vars, drop = FALSE]
          }

          # Stop if no columns remain after filtering
          if (ncol(data) == 0) {
            shiny::showNotification("⚠️ No columns available after filtering!", type = "warning")
            return(NULL)
          }

          shiny::incProgress(0.4, detail = "Data ready for clustering...")

          # -----------------------------------------------------------
          # HCLUST
          # -----------------------------------------------------------

          # Initialize the HClust with selected options
          hc <- HClustVar$new(
            vartype = input$vartype,
            dist.metric = dist_metric,
            cah.method = input$cah_method
          )

          # Fit and cut
          shiny::incProgress(0.3, detail = "Fitting hierarchical model...")
          hc$fit(data)

          shiny::incProgress(0.3, detail = "Cutting tree...")
          hc$cut_tree(k = input$n_clusters)

          shiny::incProgress(0.2, detail = "Finalizing...")
          rv$model <- hc
          rv$model_type <- 'hclust'
          rv$clustering_done <- TRUE

          # Get the illustrative vars if any.
          if (length(illust_vars) > 0) {
            illustrative_df <-req(data_filtered())
            illustrative_df <- illustrative_df[, illust_vars, drop = FALSE]


            prediction_result <- hc$predict(illustrative_df)

            # Store the prediction as a dataframe
            rv$prediction <- data.frame(
              "Illustrative_variable" = illust_vars,
              "Variable_type" = config_df$Type[match(illust_vars, config_df$Variable)],
              "Predicted_clusters" = prediction_result$labels,
              stringsAsFactors = FALSE
            )

            # Save the matrix of similiraties.
            rv$prediction_proximities <- prediction_result$proximities

            # If no illustrative vars, set the prediction to NULL.
          } else {
            rv$prediction <- NULL
            rv$prediction_proximities <- NULL
          }

          shiny::showNotification(
            "✅ Hierarchical clustering completed successfully!",
            type = "message",
            duration = 3
          )



          # TODO: Check Kmeans
        # -----------------------------------------------------------
        # Kmeans
        # -----------------------------------------------------------

        } else if (input$algorithm == 'kmeans') {

          # Correct the typologie for kmeans
          if (input$dist_metric == "rsquare") {
            dist_metric <- "r_squared"
          } else if (input$dist_metric == "r") {
            dist_metric <- "r_signed"
          }

          # Prepare Kmeans parameters
          km_n_clusters <- input$n_clusters
          km_max_iter <- input$km_max_iter
          km_n_init <- input$km_n_init
          km_distance_metric <- dist_metric
          km_rand <- if (!is.null(input$km_random_state) && input$km_random_state > 0) input$km_random_state else NULL

          km <- KmeansVariables$new(
            n_clusters = km_n_clusters,
            max_iter = km_max_iter,
            n_init = km_n_init,
            random_state = km_rand,
            distance_metric = km_distance_metric
          )

          shiny::incProgress(0.3, detail = "Fitting k-means variables model...")
          km$fit(data)

          shiny::incProgress(0.5, detail = "Finalizing k-means results...")
          rv$model <- km
          rv$model_type <- 'kmeans'
          rv$clustering_done <- TRUE

          shiny::showNotification(
            "✅ K-means variable clustering completed successfully!",
            type = "message",
            duration = 3
          )

        } else {
          stop('Unknown algorithm selected')
        }

      }, error = function(e) {
        shiny::showNotification(
          paste("❌ Error:", e$message),
          type = "error",
          duration = 10
        )
        rv$clustering_done <- FALSE
      })
    })
  })

  # ===========================================================
  # OUTPUTS - PLOTS (HCLUST)
  # ===========================================================

  output$clustering_done <- shiny::reactive({
    return(rv$clustering_done)
  })
  shiny::outputOptions(output, "clustering_done", suspendWhenHidden = FALSE)

  # Dendrogram
  output$dendrogram_plot <- shiny::renderPlot({
    req(rv$model)
    req(rv$clustering_done)
    req(rv$model_type == 'hclust')

    tryCatch({
      if (!is.null(input$n_clusters)) {
        rv$model$plot_dendrogram(k = input$n_clusters)
      } else {
        rv$model$plot_dendrogram()
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error plotting dendrogram:", e$message), col = "red")
    })
  })

  # Aggregation levels (Elbow for HClust)
  output$agg_levels_plot <- shiny::renderPlot({
    req(rv$model)
    req(rv$clustering_done)
    req(rv$model_type == 'hclust')

    tryCatch({
      rv$model$plot_agg_levels()
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # MDS Projection (HClust)
  output$mds_plot <- shiny::renderPlot({
    req(rv$model)
    req(rv$clustering_done)
    req(rv$model_type == 'hclust')

    tryCatch({
      if (!is.null(rv$model$mds_projection)) {
        rv$model$mds_projection()
      } else if (!is.null(rv$model$plot_projection)) {
        rv$model$plot_projection()
      } else {
        plot.new()
        text(0.5, 0.5, "MDS projection not available.", cex = 1.1)
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # ===========================================================
  # OUTPUTS - PLOTS (KMEANS)
  # ===========================================================

  # Inertia Evolution Plot
  output$inertia_plot <- shiny::renderPlot({
    req(rv$model)
    req(rv$clustering_done)
    req(rv$model_type == 'kmeans')

    tryCatch({
      # Si votre modèle a une méthode plot_elbow ou inertia_evolution
      if (!is.null(rv$model$plot_elbow)) {
        rv$model$plot_elbow()
      } else if (!is.null(rv$model$inertia)) {
        # Créer un plot manuel de l'inertie
        inertia_value <- rv$model$inertia

        plot.new()
        par(mar = c(5, 5, 4, 2))
        plot(1, inertia_value,
             type = "p",
             pch = 19,
             cex = 2,
             col = "#E74C3C",
             xlab = "Iteration",
             ylab = "Inertia",
             main = "Within-Cluster Sum of Squares (Inertia)",
             cex.lab = 1.2,
             cex.main = 1.4,
             ylim = c(0, max(inertia_value) * 1.1))

        grid(col = "gray80")
        points(1, inertia_value, pch = 19, cex = 2, col = "#E74C3C")
        text(1, inertia_value,
             labels = sprintf("Inertia: %.2f", inertia_value),
             pos = 3, cex = 1.1, col = "#2C3E50")

      } else {
        plot.new()
        text(0.5, 0.5, "Inertia data not available.", cex = 1.1, col = "gray50")
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # Variable Projection (K-means)
  output$kmeans_projection_plot <- shiny::renderPlot({
    req(rv$model)
    req(rv$clustering_done)
    req(rv$model_type == 'kmeans')

    tryCatch({
      if (!is.null(rv$model$plot_projection)) {
        rv$model$plot_projection()
      } else if (!is.null(rv$model$mds_projection)) {
        rv$model$mds_projection()
      } else {
        # Créer une projection MDS manuelle si les données sont disponibles
        data <- data_filtered()
        labels <- rv$model$labels_

        if (!is.null(labels) && length(labels) == ncol(data)) {
          # Calculer la matrice de corrélation
          cor_matrix <- cor(data, use = "pairwise.complete.obs")
          dist_matrix <- as.dist(1 - abs(cor_matrix))

          # MDS
          mds_result <- cmdscale(dist_matrix, k = 2)

          # Préparer les couleurs
          colors <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12",
                      "#9B59B6", "#1ABC9C", "#E67E22", "#95A5A6")
          cluster_colors <- colors[labels]

          # Plot
          par(mar = c(5, 5, 4, 2))
          plot(mds_result[, 1], mds_result[, 2],
               col = cluster_colors,
               pch = 19,
               cex = 1.5,
               xlab = "Dimension 1",
               ylab = "Dimension 2",
               main = "MDS Projection of Variables (K-means)",
               cex.lab = 1.2,
               cex.main = 1.4)

          # Ajouter les labels
          text(mds_result[, 1], mds_result[, 2],
               labels = colnames(data),
               pos = 3,
               cex = 0.8,
               col = cluster_colors)

          # Légende
          legend("topright",
                 legend = paste("Cluster", unique(labels)),
                 col = colors[unique(labels)],
                 pch = 19,
                 cex = 0.9)

          grid(col = "gray80")
        } else {
          plot.new()
          text(0.5, 0.5, "Projection data not available.", cex = 1.1)
        }
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # Cluster Centers Heatmap
  output$centers_heatmap <- shiny::renderPlot({
    req(rv$model)
    req(rv$clustering_done)
    req(rv$model_type == 'kmeans')

    tryCatch({
      # Essayer d'obtenir les centres des clusters
      if (!is.null(rv$model$cluster_centers_)) {
        centers <- rv$model$cluster_centers_

        # Créer une heatmap des centres
        par(mar = c(5, 10, 4, 2))

        # Normaliser les centres pour la visualisation
        centers_scaled <- scale(centers)

        # Créer la heatmap
        image(t(centers_scaled),
              col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
              xlab = "Cluster",
              ylab = "",
              main = "Cluster Centers Heatmap",
              axes = FALSE,
              cex.main = 1.4)

        # Ajouter les axes
        axis(1, at = seq(0, 1, length.out = nrow(centers)),
             labels = paste("Cluster", 1:nrow(centers)),
             cex.axis = 1.1)

        axis(2, at = seq(0, 1, length.out = ncol(centers)),
             labels = colnames(centers),
             las = 2,
             cex.axis = 0.9)

        # Ajouter une légende de couleur
        legend("topright",
               legend = c("High", "Medium", "Low"),
               fill = c("#E74C3C", "white", "#3498DB"),
               cex = 0.9,
               title = "Correlation")

      } else if (!is.null(rv$model$labels_)) {
        # Alternative: montrer la composition des clusters
        data <- data_filtered()
        labels <- rv$model$labels_

        # Calculer les centres manuellement
        centers <- matrix(NA, nrow = input$n_clusters, ncol = ncol(data))
        colnames(centers) <- colnames(data)
        rownames(centers) <- paste("Cluster", 1:input$n_clusters)

        for (k in 1:input$n_clusters) {
          vars_in_cluster <- which(labels == k)
          if (length(vars_in_cluster) > 0) {
            cluster_data <- data[, vars_in_cluster, drop = FALSE]
            # Calculer la corrélation moyenne avec toutes les variables
            for (j in 1:ncol(data)) {
              cors <- cor(data[, j], cluster_data, use = "pairwise.complete.obs")
              centers[k, j] <- mean(abs(cors), na.rm = TRUE)
            }
          }
        }

        # Plot heatmap
        par(mar = c(5, 10, 4, 2))
        image(t(centers),
              col = colorRampPalette(c("white", "#3498DB", "#E74C3C"))(100),
              xlab = "Cluster",
              ylab = "",
              main = "Cluster-Variable Correlation Strength",
              axes = FALSE,
              cex.main = 1.4)

        axis(1, at = seq(0, 1, length.out = nrow(centers)),
             labels = rownames(centers),
             cex.axis = 1.1)

        axis(2, at = seq(0, 1, length.out = ncol(centers)),
             labels = colnames(centers),
             las = 2,
             cex.axis = 0.9)

      } else {
        plot.new()
        text(0.5, 0.5, "Cluster centers not available.", cex = 1.1)
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # ===========================================================
  # OUTPUTS - COMMON PLOTS
  # ===========================================================

  # Silhouette (commun aux deux algorithmes)
  output$silhouette_plot <- shiny::renderPlot({
    req(rv$model)
    req(rv$clustering_done)

    tryCatch({
      # Try plotting silhouette (if available) for either model
      if (!is.null(rv$model$plot_silhouette)) {
        rv$model$plot_silhouette()
      } else if (!is.null(rv$model$silhouette)) {
        # if a silhouette() returns data, try to plot basic silhouette
        sil <- rv$model$silhouette()
        if (!is.null(sil)) {
          plot(sil)
        } else {
          plot.new()
          text(0.5, 0.5, "No silhouette available.")
        }
      } else {
        plot.new()
        text(0.5, 0.5, "Silhouette plot not available for selected algorithm.", cex = 1.1)
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # ===========================================================
  # OUTPUTS - TABLES
  # ===========================================================

  # -----------------------------------------------------------
  # OUTPUTS - Prediction
  # -----------------------------------------------------------

  output$prediction_result <- DT::renderDT({
    req(rv$prediction)

    prediction <- rv$prediction

    DT::datatable(
      prediction,
      options = list(
        pageLength = 20,
        dom = 't'
      ),
      rownames = FALSE,
      class = 'cell-border stripe'
    ) %>%
      DT::formatStyle(
        'Predicted_clusters',
        backgroundColor = '#E8F5E9',
        fontWeight = 'bold'
      )
  })

  # Prediction proximities table
  output$prediction_proximities_table <- DT::renderDT({
    req(rv$prediction_proximities)

    proximities_df <- rv$prediction_proximities

    DT::datatable(
      proximities_df,
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      class = 'cell-border stripe compact'
    ) %>%
      # Color code based on proximity values
      DT::formatStyle(
        columns = 2:ncol(proximities_df),
        backgroundColor = DT::styleInterval(
          cuts = c(0.3, 0.5, 0.7, 0.85),
          values = c('#FFEBEE', '#FFF3E0', '#FFF9C4', '#E8F5E9', '#C8E6C9')
        )
      ) %>%
      # Format numbers with 3 decimals
      DT::formatRound(columns = 2:ncol(proximities_df), digits = 3) %>%
      # Bold for high proximity (> 0.7)
      DT::formatStyle(
        columns = 2:ncol(proximities_df),
        fontWeight = DT::styleInterval(
          cuts = c(0.7),
          values = c('normal', 'bold')
        )
      )
  })

  # -----------------------------------------------------------
  # OUTPUTS - Summary tables
  # -----------------------------------------------------------
  output$cluster_summary_table <- DT::renderDT({
    req(rv$model)
    req(rv$clustering_done)

    summary_data <- tryCatch(rv$model$summary(), error = function(e) NULL)
    if (is.null(summary_data)) {
      DT::datatable(data.frame())
    } else {
      # If clust_summary exists, show it. Else try to show whole summary object
      df <- if (!is.null(summary_data$clust_summary)) summary_data$clust_summary else as.data.frame(summary_data)
      dt <- DT::datatable(df, options = list(pageLength = 20, dom = 't'), rownames = FALSE)
      dt
    }
  })

  output$cluster_members_table <- DT::renderDT({
    req(rv$model)
    req(rv$clustering_done)

    summary_data <- tryCatch(rv$model$summary(), error = function(e) NULL)
    if (is.null(summary_data)) {
      DT::datatable(data.frame())
    } else {
      df <- if (!is.null(summary_data$clust_members)) summary_data$clust_members else as.data.frame(summary_data)
      dt <- DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE), rownames = TRUE)
      if ('own_cluster_R2' %in% colnames(df)) {
        dt <- dt %>% DT::formatStyle('own_cluster_R2', backgroundColor = DT::styleInterval(c(0.5, 0.7), c('#FADBD8', '#FCF3CF', '#D5F4E6')))
      }
      dt
    }
  })

  # Average Silhouette Display
  output$avg_silhouette_display <- shiny::renderUI({
    req(rv$model)
    req(rv$clustering_done)

    summary_data <- tryCatch(rv$model$summary(), error = function(e) NULL)
    avg_sil <- summary_data$avg_silhouette

    if (is.null(avg_sil)) {
      shiny::tags$div(
        style = "color: #95A5A6; font-size: 14px;",
        "Not available"
      )
    } else {
      # Color based on silhouette value
      color <- if (avg_sil >= 0.7) {
        "#4CAF50"  # Green - excellent
      } else if (avg_sil >= 0.5) {
        "#FFC107"  # Yellow - good
      } else if (avg_sil >= 0.25) {
        "#FF9800"  # Orange - fair
      } else {
        "#F44336"  # Red - poor
      }

      quality <- if (avg_sil >= 0.7) {
        "Excellent"
      } else if (avg_sil >= 0.5) {
        "Good"
      } else if (avg_sil >= 0.25) {
        "Fair"
      } else {
        "Poor"
      }

      shiny::tagList(
        shiny::tags$div(
          style = paste0("font-size: 48px; font-weight: bold; color: ", color, ";"),
          sprintf("%.3f", avg_sil)
        ),
        shiny::tags$div(
          style = paste0("font-size: 16px; color: ", color, "; margin-top: 5px;"),
          quality
        )
      )
    }
  })

  # Total Variance Explained Display
  output$total_var_explained_display <- shiny::renderUI({
    req(rv$model)
    req(rv$clustering_done)

    summary_data <- tryCatch(rv$model$summary(), error = function(e) NULL)
    total_var <- summary_data$total_var_explained

    if (is.null(total_var)) {
      shiny::tags$div(
        style = "color: #95A5A6; font-size: 14px;",
        "Not available"
      )
    } else {
      # Convert to percentage if needed
      if (total_var <= 1) {
        total_var <- total_var * 100
      }

      # Color based on variance explained
      color <- if (total_var >= 70) {
        "#4CAF50"  # Green - excellent
      } else if (total_var >= 50) {
        "#FFC107"  # Yellow - good
      } else if (total_var >= 30) {
        "#FF9800"  # Orange - fair
      } else {
        "#F44336"  # Red - poor
      }

      quality <- if (total_var >= 70) {
        "Excellent"
      } else if (total_var >= 50) {
        "Good"
      } else if (total_var >= 30) {
        "Fair"
      } else {
        "Poor"
      }

      shiny::tagList(
        shiny::tags$div(
          style = paste0("font-size: 48px; font-weight: bold; color: ", color, ";"),
          sprintf("%.1f%%", total_var)
        ),
        shiny::tags$div(
          style = paste0("font-size: 16px; color: ", color, "; margin-top: 5px;"),
          quality
        )
      )
    }
  })

  # Centroids Correlation Table
  output$centroids_correlation_table <- DT::renderDT({
    req(rv$model)
    req(rv$clustering_done)

    summary_data <- tryCatch(rv$model$summary(), error = function(e) NULL)
    centroids_df <- summary_data$centroids_correlations

    if (is.null(centroids_df)) {
      # Return empty datatable with message
      DT::datatable(
        data.frame(Message = "Centroids correlation matrix not available"),
        options = list(dom = 't'),
        rownames = FALSE
      )
    } else {

      DT::datatable(
        centroids_df,
        options = list(
          pageLength = 10,
          dom = 't',
          scrollX = TRUE
        ),
        rownames = FALSE
      ) %>%
        # Format correlation values
        DT::formatRound(columns = 3:ncol(centroids_df), digits = 3)
    }
  })

  # All Clusters R² Table
  output$all_clusters_R2_table <- DT::renderDT({
    req(rv$model)
    req(rv$clustering_done)

    summary_data <- tryCatch(rv$model$summary(), error = function(e) NULL)

    if (is.null(summary_data) || is.null(summary_data$all_clusters_R2)) {
      # Return empty datatable with message
      DT::datatable(
        data.frame(Message = "R² matrix not available"),
        options = list(dom = 't'),
        rownames = FALSE
      )
    } else {
      all_clusters_R2 <- summary_data$all_clusters_R2

      DT::datatable(
        all_clusters_R2,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          scrollY = "500px",
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        rownames = FALSE,
        class = 'cell-border stripe compact'
      ) %>%
        # Colorer selon les valeurs de R²
        DT::formatStyle(
          columns = 2:ncol(all_clusters_R2),
          backgroundColor = DT::styleInterval(
            cuts = c(0.3, 0.5, 0.7, 0.85),
            values = c('#FFEBEE', '#FFF3E0', '#FFF9C4', '#E8F5E9', '#C8E6C9')
          )
        ) %>%
        # Formater les nombres avec 3 décimales
        DT::formatRound(columns = 2:ncol(all_clusters_R2), digits = 3) %>%
        # Mettre en gras les valeurs > 0.7 (forte corrélation)
        DT::formatStyle(
          columns = 2:ncol(all_clusters_R2),
          fontWeight = DT::styleInterval(
            cuts = c(0.7),
            values = c('normal', 'bold')
          )
        )
    }
  })


  # ===========================================================
  # DOWNLOAD HANDLER
  # ===========================================================

  output$download_results <- shiny::downloadHandler(
    filename = function() {
      if (!is.null(rv$model_type) && rv$model_type == 'kmeans') {
        paste0("kmeansvar_results_", Sys.Date(), ".csv")
      } else {
        paste0("hclustvar_results_", Sys.Date(), ".csv")
      }
    },
    content = function(file) {
      req(rv$model)
      req(rv$clustering_done)

      summary_data <- tryCatch(rv$model$summary(), error = function(e) NULL)
      if (!is.null(summary_data) && !is.null(summary_data$clust_members)) {
        write.csv(summary_data$clust_members, file, row.names = TRUE)
      } else if (!is.null(summary_data)) {
        # fallback: write the whole summary as csv if possible
        write.csv(as.data.frame(summary_data), file, row.names = TRUE)
      } else {
        write.csv(data.frame(), file)
      }
    }
  )
}
