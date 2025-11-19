#' Launch HClustVar Shiny Application
#'
#' @description
#' Launches an interactive Shiny application for hierarchical variable clustering.
#'
#' @param ... Additional arguments passed to \code{shiny::runApp()}
#'
#' @return No return value, launches Shiny app
#'
#' @examples
#' \dontrun{
#' # Launch the app
#' run_hclustvar_app()
#' }
#'
#' @export
varclust_gui <- function(...) {

  # Vérifier que shiny est installé
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' required. Install it with: install.packages('shiny')")
  }

  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' required. Install it with: install.packages('DT')")
  }

  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' required. Install it with: install.packages('readxl')")
  }

  # Créer l'application
  app <- shiny::shinyApp(
    ui = app_ui(),
    server = app_server
  )

  # Lancer l'application
  shiny::runApp(app, ...)
}
