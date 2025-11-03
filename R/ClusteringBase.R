#' Base class for variable clustering algorithms
#'
#' @description
#' Defines the common interface and core attributes used by all clustering algorithm classes.
#' This class is **abstract** and must not be instantiated directly.
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
ClusteringBase <- R6::R6Class(

  # Calling string.
  "ClusteringBase",


  # ==========================================================================
  # PRIVATE FIELDS
  # ==========================================================================
  private = list(

    .data = NULL,
    .n_clusters = NULL,
    .labels = NULL,
    .fitted = FALSE,

    # Indices of variable types (computed once)
    .quanti_indices = NULL,
    .quali_indices = NULL,

    # -----------------------------------------------------------------------
    # CORE VALIDATION METHODS
    # -----------------------------------------------------------------------

    #' Validate input data structure
    #' @param data Data.frame or matrix
    #' @return The data as a dataframe and column named if needed.
    validate_data_structure = function(data) {

      # Check type
      if (!is.data.frame(data) && !is.matrix(data)) {
        stop("'data' must be a data.frame or matrix, got: ", class(data)[1])
      }

      # Convert matrix to data.frame for consistency
      if (is.matrix(data)) {
        data <- as.data.frame(data)
      }

      # Check dimensions
      if (ncol(data) < 2) {
        stop("At least 2 variables are required for clustering, got: ", ncol(data))
      }

      if (nrow(data) < 2) {
        stop("At least 2 observations are required, got: ", nrow(data))
      }

      # Check for column names
      if (is.null(colnames(data))) {
        warning("No column names provided. Using default names: V1, V2, ...")
        colnames(data) <- paste0("V", 1:ncol(data))
      }
      return(data)
    },


    #' Validate number of clusters
    #' @param n_clusters Integer or NULL
    validate_n_clusters = function(n_clusters) {

      # Return NULL if needed as some algorithms dosn't require a specified k.
      if (is.null(n_clusters)) {
        return(NULL)
      }

      # Check type and single element numeric vector.
      if (!is.numeric(n_clusters) || length(n_clusters) != 1) {
        stop("'n_clusters' must be a single integer, got: ",
             paste(n_clusters, collapse = ", "))
      }

      # Convert into integer if float.
      n_clusters <- as.integer(n_clusters)

      if (n_clusters < 1) {
        stop("'n_clusters' must be at least 1, got: ", n_clusters)
      }

      if (!is.null(private$.data) && n_clusters > ncol(private$.data)) {
        stop("'n_clusters' (", n_clusters,
             ") cannot exceed number of variables (", ncol(private$.data), ")")
      }

      return(n_clusters)
    },

    #' Check for missing values. Raise an error if found one.
    check_missing_values = function() {

      # Count missing values per variable
      missing_counts <- colSums(is.na(private$.data))

      # Check if there any NA in the dataframe.
      total_missing <- sum(missing_counts)

      if (total_missing > 0) {
        stop("Your data contains NA values: ", missing_counts)
      }
    },


    #' Detect variable types (quantitative vs qualitative)
    detect_variable_types = function() {

      # Get the quantitative and the qualitatives data indices.
      quant.index <- unname(which(sapply(private$.data, is.numeric)))
      qual.index <- setdiff(seq_along(private$.data), quant.index)

      return(list(quant.index, qual.index))
    }
  ),  # end private


  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list (

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------

    initialize = function() {
      # Constructor preventing direct implementation -> Abstract class.
      if (class(self)[1] == "ClusteringBase") {
        stop("ClusteringBase is an abstract class and shouldn't be directly instanciated.")
      }
    },

    # -----------------------------------------------------------------------
    # Abstract methods (to be implemented by child classes)
    # -----------------------------------------------------------------------

    fit = function() {
      stop("The fit() method should be implemented in a child class.")
    },

    predict = function() {
      stop("The predict() method should be implemented in a child class.")
    },

    summary = function() {
      stop("The summary() method should be implemented in a child class.")
    },

    print = function() {
      stop("The print() method should be implemented in a child class.")
    },

    # -----------------------------------------------------------------------
    # Utility methods (available to all child classes)
    # -----------------------------------------------------------------------


    #' Load and validate clustering input data
    #'
    #' @description
    #' Performs all necessary checks and preprocessing steps to ensure the input
    #' dataset is compatible with the clustering algorithm.
    #'
    #' This method is intended to be called inside child classes (e.g., within
    #' their `fit()` method) before running the algorithm. It centralizes the
    #' validation logic: verifying data structure, missing values, and detecting
    #' variable types (quantitative vs qualitative).
    load_and_check_data = function(data) {

      if (!is.null(data)) {
        # Validate and save data
        private$.data <- private$validate_data_structure(data)
        private$check_missing_values()

        types <- private$detect_variable_types()
        private$.quanti_indices <- types[[1]]
        private$.quali_indices <- types[[2]]

      }
    },

    #' Get quantitative variables only
    #' @return Data.frame with quantitative variables
    get_quanti_data = function() {
      if (length(private$.quanti_indices) == 0) {
        stop("No quantitative variables found in data.")
      }
      return(private$.data[, private$.quanti_indices, drop = FALSE]) # Drop to keep the dataframe format.
    },


    #' Get qualitative variables only
    #' @return Data.frame with qualitative variables
    get_quali_data = function() {
      if (length(private$.quali_indices) == 0) {
        stop("No qualitative variables found in data.")
      }
      return(private$.data[, private$.quali_indices, drop = FALSE]) # Drop to keep the dataframe format.
    },

    #' Validate data compatibility for specific algorithms.
    #' To be be used within child classes to validate data types.
    #' @param type Character: "quant", "qual", "mixed".
    validate_algorithm_requirements = function(type) {

      has_quanti <- length(private$.quanti_indices) > 1
      has_quali <- length(private$.quali_indices) > 1


      if (!type %in% c("quant", "qual", "mixed")) {
        stop("Invalid requirement type.")
      }

      if (type == "quant" && !has_quanti) {
        stop("This algorithm requires at least one quantitative variable.")
      }

      if (type == "qual" && !has_quali) {
        stop("This algorithm requires at least one qualitative variable.")
      }
    }

  ), # End public

  # ==========================================================================
  # ACTIVE BINDINGS (getters/setters)
  # ==========================================================================

  active = list(

    # -------------------------------------------------------------------------
    # data getter/setter
    # -------------------------------------------------------------------------
    data = function(value) {

      if (missing(value)) {
        return(private$.data)
      } else {
        private$.data <- value
      }
    },

    # -------------------------------------------------------------------------
    # n_clusters getter/setter
    # -------------------------------------------------------------------------
    n_clusters = function(value) {
      if (missing(value)) {
        return(private$.n_clusters)
      }

      # Check that the value is a correct integer.
      if (!is.numeric(value) || length(value) != 1 || value <= 0 || value %% 1 != 0) {
        stop("'n_clusters' must be a positive integer scalar.")
      }

      private$.n_clusters <- private$validate_n_clusters(value)
    },

    # -------------------------------------------------------------------------
    # labels getter/setter
    # -------------------------------------------------------------------------
    labels = function(value) {
      if (missing(value)) {
        return(private$.labels)
      }

      # Cheking length consistency with data.
      if (!is.null(private$.data) && length(value) != nrow(private$.data)) {
        stop("'labels' length must match number of observations.")
      }

      private$.labels <- value
    },

    # -------------------------------------------------------------------------
    # fitted getter/setter
    # -------------------------------------------------------------------------
    fitted = function(value) {
      if (missing(value)) {
        return(private$.fitted)
      }

      # Must be boolean.
      if (!is.logical(value) || length(value) != 1) {
        stop("'fitted' must be a single logical value (TRUE/FALSE).")
      }

      private$.fitted <- value
    },

    # -------------------------------------------------------------------------
    # Quantitative / qualitative indices getter
    # -------------------------------------------------------------------------
    quanti_indices = function() private$.quanti_indices,
    quali_indices  = function() private$.quali_indices
  )
)



