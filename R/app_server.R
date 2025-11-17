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
    hclust_model = NULL,
    clustering_done = FALSE
  )

  # ===========================================================
  # DATA LOADING
  # ===========================================================

  # Charger les données depuis le fichier uploadé
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

        # Convertir en data.frame
        data <- as.data.frame(data)

        # Détecter les colonnes qualitatives et les convertir en factor
        for (col in names(data)) {
          if (is.character(data[[col]]) ||
              (is.numeric(data[[col]]) && length(unique(data[[col]])) <= 10)) {
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

  output$data_info <- shiny::renderUI({
    req(data_loaded())

    data <- data_loaded()
    n_rows <- nrow(data)
    n_cols <- ncol(data)
    n_numeric <- sum(sapply(data, is.numeric))
    n_factor <- sum(sapply(data, is.factor))

    shiny::wellPanel(
      style = "background-color: #D5F4E6; border-left: 5px solid #27AE60;",
      shiny::h4("✅ Data loaded successfully!", style = "color: #27AE60; margin-top: 0;"),
      shiny::tags$ul(
        shiny::tags$li(paste("Rows:", n_rows)),
        shiny::tags$li(paste("Columns:", n_cols)),
        shiny::tags$li(paste("Numeric variables:", n_numeric)),
        shiny::tags$li(paste("Categorical variables:", n_factor))
      )
    )
  })

  output$data_table <- DT::renderDT({
    req(data_loaded())

    DT::datatable(
      data_loaded(),
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

    req(data_loaded())

    shiny::withProgress(message = 'Running clustering...', value = 0, {

      tryCatch({

        # Étape 1: Initialisation
        shiny::incProgress(0.2, detail = "Initializing model...")

        dist_metric <- if (input$vartype == "quant") {
          input$dist_metric
        } else {
          NULL
        }

        hc <- HClustVar$new(
          vartype = input$vartype,
          dist.metric = dist_metric,
          cah.method = input$cah_method
        )

        # Étape 2: Fit
        shiny::incProgress(0.3, detail = "Fitting model...")
        hc$fit(data_loaded())

        # Étape 3: Cut tree
        shiny::incProgress(0.3, detail = "Cutting tree...")
        hc$cut_tree(k = input$n_clusters)

        # Étape 4: Stockage
        shiny::incProgress(0.2, detail = "Finalizing...")
        rv$hclust_model <- hc
        rv$clustering_done <- TRUE

        shiny::showNotification(
          "✅ Clustering completed successfully!",
          type = "message",
          duration = 3
        )

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
  # OUTPUTS - PLOTS
  # ===========================================================

  output$clustering_done <- shiny::reactive({
    return(rv$clustering_done)
  })
  shiny::outputOptions(output, "clustering_done", suspendWhenHidden = FALSE)

  # Dendrogram
  output$dendrogram_plot <- shiny::renderPlot({
    req(rv$hclust_model)
    req(rv$hclust_model$fitted)

    rv$hclust_model$plot_dendrogram(k = input$n_clusters)
  })

  # Aggregation levels
  output$agg_levels_plot <- shiny::renderPlot({
    req(rv$hclust_model)
    req(rv$hclust_model$fitted)

    rv$hclust_model$plot_agg_levels()
  })

  # Silhouette
  output$silhouette_plot <- shiny::renderPlot({
    req(rv$hclust_model)
    req(rv$clustering_done)

    tryCatch({
      rv$hclust_model$plot_silhouette()
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # MDS
  output$mds_plot <- shiny::renderPlot({
    req(rv$hclust_model)
    req(rv$clustering_done)

    tryCatch({
      rv$hclust_model$mds_projection()
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), col = "red")
    })
  })

  # ===========================================================
  # OUTPUTS - TABLES
  # ===========================================================

  output$cluster_summary_table <- DT::renderDT({
    req(rv$hclust_model)
    req(rv$clustering_done)

    summary_data <- rv$hclust_model$summary()

    DT::datatable(
      summary_data$clust_summary,
      options = list(
        pageLength = 20,
        dom = 't'
      ),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        'prop_explained',
        backgroundColor = DT::styleInterval(
          c(0.5, 0.7),
          c('#FADBD8', '#FCF3CF', '#D5F4E6')
        )
      )
  })

  output$cluster_members_table <- DT::renderDT({
    req(rv$hclust_model)
    req(rv$clustering_done)

    summary_data <- rv$hclust_model$summary()

    DT::datatable(
      summary_data$clust_members,
      options = list(
        pageLength = 20,
        scrollX = TRUE
      ),
      rownames = TRUE
    ) %>%
      DT::formatStyle(
        'own_cluster_R2',
        backgroundColor = DT::styleInterval(
          c(0.5, 0.7),
          c('#FADBD8', '#FCF3CF', '#D5F4E6')
        )
      )
  })

  # ===========================================================
  # DOWNLOAD HANDLER
  # ===========================================================

  output$download_results <- shiny::downloadHandler(
    filename = function() {
      paste0("hclustvar_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$hclust_model)
      req(rv$clustering_done)

      summary_data <- rv$hclust_model$summary()
      write.csv(summary_data$clust_members, file, row.names = TRUE)
    }
  )
}
