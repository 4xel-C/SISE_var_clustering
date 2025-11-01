HClustVar <- R6::R6Class(
  "HClustVar",
  inherit = ClusteringBase,


  # ==========================================================================
  # PRIVATE FIELDS
  # ==========================================================================
  private = list(
    .dendrogram = NULL,
    .height_cut = NULL,
    .dissimilarity_matrix = NULL,

    # -----------------------------------------------------------------------
    # CORE VALIDATION METHODS
    # -----------------------------------------------------------------------
    #' Validate height_cut input.
    #' @param height Integer or NULL
    validate_height_cut = function(height) {

      # Return NULL if no value imputed.
      if (is.null(height)) {
        return(NULL)
      }

      # Check type and single element numeric vector.
      if (!is.numeric(height) || length(height) != 1) {
        stop("'height' must be a single integer, got: ",
             paste(height, collapse = ", "))
      }

      # Convert into integer if float.
      height <- as.integer(height)

      if (height < 1) {
        stop("'height' must be at least 1, got: ", height)
      }
      return(height)
    }
  ),

  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list(

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------
    initialize = function(data, n_clusters = NULL, height_cut = NULL) {
      super$initialize(data, n_clusters)

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
