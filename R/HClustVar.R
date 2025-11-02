
#' Variable Clustering class based on AHC with distances matrix.
#'
#' @description
#'
#'
#' @details
#' The `ClusteringBase` class provides:
#' - private fields for storing data, number of clusters, labels, and fit status
#' - Input validation and data preprocessing utilities
#' - Detection of variable types (quantitative vs qualitative)
#' - Abstract methods that must be implemented by child classes
#'
#' @note
#' This class is not exported outside the package. It serves only as an internal prototype.
#'
#' @noRd
HClustVar <- R6::R6Class(
  "HClustVar",
  inherit = ClusteringBase,


  # ==========================================================================
  # PRIVATE FIELDS
  # ==========================================================================
  private = list(

    # Metric to calculate distances between variables.
    .dist.metric = NULL

  ),

  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list(

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------
    initialize = function(dist.metric="auto") {

      # Checking metric select.
      if (dist.metric %in% c("rsquare", "r", "disimilarity", "auto")) {
        stop("Parameter 'dist.metric' has invalid value. Choose: 'rsquare', 'r', 'disimilarity', 'auto'. Got: ", dist.metric)
      }

      private$.dist.metric <- dist.metric
    },

    fit = function() {
      # TODO: A implémenter
    },

    cut_tree = function(k = NULL, h = NULL) {
      # TODO: A implémenter
    },

    # Dendrogramme
    plot_dendrogram = function(k = NULL, ...) {
      # TODO: A implémenter
    }
  ),

  # ==========================================================================
  # Getters/Setters
  # ==========================================================================
  active = list(
    dissimilarity = function() {private$.dissimilarity_matrix}
  )
)
