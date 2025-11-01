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
      quant.index <- which(sapply(private$.data, is.numeric))
      qual.index <- setdiff(seq_along(private$.data), quant.index)

      return(list(quant.index, qual.index))
    },


    #' Validate data compatibility for specific algorithms.
    #' To be be used while instantiating child classes to validate data types.
    #' @param requires_type Character: "quanti", "quali", or "mixed"
    validate_algorithm_requirements = function(requires_type = "mixed") {

      has_quanti <- length(private$.quanti_indices) > 0
      has_quali <- length(private$.quali_indices) > 0

      if (requires_type == "quanti" && !has_quanti) {
        stop("This algorithm requires at least one quantitative variable.")
      }

      if (requires_type == "quali" && !has_quali) {
        stop("This algorithm requires at least one qualitative variable.")
      }

      if (requires_type == "mixed" && !(has_quanti && has_quali)) {
        warning("Algorithm expects mixed data, but data is homogeneous (",
                ifelse(has_quanti, "all quantitative", "all qualitative"), ")")
      }
    }
  ),  # end private


  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list (

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------

    initialize = function(data = NULL, n_clusters = NULL) {

      # Constructor preventing direct implementation -> Abstract class.
      if (class(self)[1] == "ClusteringBase") {
        stop("ClusteringBase is an abstract class and shouldn't be directly instanciated.")
      }

      if (!is.null(data)) {
        # Validate and save data
        private$.data <- private$validate_data_structure(data)
        private$.n_clusters <- n_clusters
        private$check_missing_values()

        types <- private$detect_variable_types()
        private$.quanti_indices <- types[[1]]
        private$.quali_indices <- types[[2]]
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
    }
  ), # End public

  # ==========================================================================
  # ACTIVE BINDINGS (getters/setters)
  # ==========================================================================

  active = list(

    # Data getter (read-only from outside)
    data = function() {return(private$.data)},

    # n_clusters getter
    n_clusters = function(value) {return(private$.n_clusters)},

    # labels getter (read-only)
    labels = function() {return(private$.labels)},

    # Quantitative variable indices
    quanti_indices = function() {return(private$.quanti_indices)},

    # Qualitative variable indices
    quali_indices = function() {return(private$.quali_indices)}

  ) # end active
)


