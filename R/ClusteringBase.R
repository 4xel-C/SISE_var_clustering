#' @title Abstract Base Class for Variable Clustering Algorithms
#'
#' @description
#' The `ClusteringBase` class defines the **foundation for variable clustering algorithms**.
#' It provides input validation, preprocessing utilities, and type detection logic
#' (quantitative vs qualitative variables).
#'
#' This class is **abstract** and should **not be directly instantiated**.
#' It is meant to be inherited by concrete classes implementing specific
#' clustering algorithms.
#'
#' @details
#' The `ClusteringBase` R6 class offers:
#'
#' - **Private fields** for storing:
#'   - `.data` : the input dataset.
#'   - `.n_clusters` : the number of clusters.
#'   - `.labels` : vector of cluster assignments.
#'   - `.fitted` : logical flag indicating fit status.
#'   - `.quanti_indices` / `.quali_indices` : detected quantitative and qualitative variable indices.
#'
#' - **Core validation private methods**:
#'   - Input structure checking (`validate_data_structure()`): Validate dataframe/matrix type, and check for number of lines/columns.
#'   - Cluster count validation (`validate_n_clusters()`): Validate the number of clusters for setter.
#'   - Missing value detection (`check_missing_values()`): Check for missing values in the dataframe and raise an error if found.
#'   - Variable type detection (`detect_variable_types()`): Check for variable type and save the indices in the class (.quanti_indices, .quali_indices).
#'
#' - **Public utilities**:
#'   - Data loading and checking (`load_and_check_data()`): Load the dataframe in the class and run all the checking methods.
#'   - Extracting subsets (`get_quanti_data()`, `get_quali_data()`): return a subset of the dataframe containing all quantitative/qualitative variables.
#'   - Compatibility checking for algorithm types (`validate_algorithm_requirements()`): Compare the datatype on which to execute the algorithm with the loaded dataframe.
#'
#' - **Abstract methods** to be overridden:
#'   - `fit()`
#'   - `predict()`
#'   - `summary()`
#'   - `print()`
#'
#' Attempting to instantiate this class directly will raise an error.
#'
#' @note
#' This is an **internal abstract base class**.
#' It should not be exported or directly accessible to package users.
#'
#' @examples
#' \dontrun{
#' # Example: inheriting from ClusteringBase
#' MyAlgo <- R6::R6Class(
#'   "MyAlgo",
#'   inherit = ClusteringBase,
#'   public = list(
#'     fit = function(data) {
#'       self$load_and_check_data(data)
#'       message("Algorithm fitting logic here.")
#'     }
#'   )
#' )
#'
#' obj <- MyAlgo$new()
#' obj$fit(iris[, 1:4])
#' }
#'
#' @keywords internal
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
    #' @return The data as a dataframe with column named if needed.
    #' @keywords internal
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
    #' @keywords internal
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

    #' Check for missing values in the data. Raise an error if found one.
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
    #' @return The list containing vector of quantitative indices, and qualitative indices.
    #' @keywords internal
    detect_variable_types = function() {

      is_quanti <- sapply(private$.data, function(x) {
        is.numeric(x) && !is.factor(x)
      })

      quant_index <- which(is_quanti)
      qual_index <- setdiff(seq_along(private$.data), quant_index)

      return(list(
        quanti = unname(quant_index),
        quali = unname(qual_index)
      ))
    }
  ),  # end private


  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list (

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------

    #' @description
    #' Constructor for the `ClusteringBase` class.
    #'
    #' This constructor prevents direct instantiation of the abstract base class.
    #' Child classes inheriting from `ClusteringBase` can implement
    #' their own constructor if needed.
    #'
    #' @examples
    #' \dontrun{
    #' MyAlgo <- R6::R6Class(
    #'   "MyAlgo",
    #'   inherit = ClusteringBase,
    #'   public = list(
    #'     initialize = function() {
    #'       message("MyAlgo initialized.")
    #'     }
    #'   )
    #' )
    #' }
    #' @keywords internal
    initialize = function() {

      # Prevent direct instantiation.
      if (class(self)[1] == "ClusteringBase") {
        stop("ClusteringBase is an abstract class and shouldn't be directly instanciated.")
      }
    },

    # -----------------------------------------------------------------------
    # Abstract methods (to be implemented by child classes)
    # -----------------------------------------------------------------------

    #' @description
    #' Abstract method to fit the clustering algorithm on the data (dataframe or matrix).
    #'
    #' @details
    #' This method must be implemented by any subclass inheriting from
    #' `ClusteringBase`. It defines the core algorithm for fitting the clustering
    #' model on the input data.
    #'
    #' @noRd
    fit = function() {
      stop("The fit() method should be implemented in a child class.")
    },

    #' @description
    #' Abstract method to predict cluster assignments for new data.
    #'
    #' @details
    #' To be implemented by subclasses. Defines how to assign clusters to new
    #' observations based on the fitted model.
    #'
    #' @noRd
    predict = function() {
      stop("The predict() method should be implemented in a child class.")
    },

    #' @description
    #' Abstract summary method.
    #'
    #' @details
    #' Should display a human-readable summary of the clustering model
    #' (e.g., number of clusters, inertia, etc.).
    #'
    #' @noRd
    summary = function() {
      stop("The summary() method should be implemented in a child class.")
    },

    #' @description
    #' Abstract print method.
    #'
    #' @details
    #' To be implemented by subclasses to control how objects of the class
    #' are printed or displayed.
    #'
    #' @noRd
    print = function() {
      stop("The print() method should be implemented in a child class.")
    },

    # -----------------------------------------------------------------------
    # Utility methods (available to all child classes)
    # -----------------------------------------------------------------------


    #' @title Load and validate clustering input data
    #'
    #' @description
    #' Performs all necessary checks and preprocessing steps to ensure that
    #' the input dataset is compatible with the clustering algorithm.
    #'
    #' This method is intended to be called **inside child classes** (for example,
    #' within their `fit()` method) before executing the algorithm.
    #' It centralizes all validation logic:
    #'
    #' - verifies data structure and dimensions
    #' - checks for missing values
    #' - detects variable types (quantitative vs qualitative)
    #'
    #' @param data A `data.frame` or `matrix` containing the variables to cluster.
    #'
    #' @return
    #' Invisibly updates the internal `.data`, `.quanti_indices`, and `.quali_indices`
    #' fields of the object.
    #'
    #' @examples
    #' \dontrun{
    #' my_algo <- MyAlgo$new()
    #' my_algo$load_and_check_data(iris[, 1:4])
    #' }
    #'
    #' @noRd
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

    #' @title Get quantitative variables only
    #'
    #' @description
    #' Returns the subset of the dataset containing **quantitative variables** only.
    #'
    #' This method can be used by subclasses (or internally) to access the numeric
    #' variables detected during the `load_and_check_data()` validation process.
    #'
    #' @return
    #' A `data.frame` containing only the quantitative (numeric) variables.
    #'
    #' @examples
    #' \dontrun{
    #' my_algo <- MyAlgo$new()
    #' my_algo$load_and_check_data(iris)
    #' quanti_data <- my_algo$get_quanti_data()
    #' head(quanti_data)
    #' }
    #'
    #' @keywords internal
    get_quanti_data = function() {
      if (length(private$.quanti_indices) == 0) {
        stop("No quantitative variables found in data.")
      }
      return(private$.data[, private$.quanti_indices, drop = FALSE]) # Drop to keep the dataframe format.
    },



    #' @title Get qualitative variables only
    #'
    #' @description
    #' Returns the subset of the dataset containing **qualitative variables** only.
    #'
    #' This method can be used by subclasses (or internally) to access the categorical
    #' variables detected during the `load_and_check_data()` validation process.
    #'
    #' @return
    #' A `data.frame` containing only the qualitative (categorical or factor) variables.
    #'
    #' @examples
    #' \dontrun{
    #' my_algo <- MyAlgo$new()
    #' my_algo$load_and_check_data(iris)
    #' quali_data <- my_algo$get_quali_data()
    #' head(quali_data)
    #' }
    #'
    #' @noRd
    get_quali_data = function() {
      if (length(private$.quali_indices) == 0) {
        stop("No qualitative variables found in data.")
      }
      return(private$.data[, private$.quali_indices, drop = FALSE]) # Drop to keep the dataframe format.
    },



    #' @title Validate data compatibility for specific clustering algorithms
    #'
    #' @description
    #' Checks that the loaded dataset meets the requirements of the clustering
    #' algorithm regarding the type and number of variables.
    #'
    #' This method is designed to be called from child classes before running
    #' their algorithm (usually inside the `fit()` method), ensuring that the
    #' expected data types (quantitative, qualitative, or mixed) are available.
    #'
    #' @param type Character string specifying the expected variable type:
    #'   - `"quant"` for algorithms requiring at least two quantitative variables.
    #'   - `"qual"` for algorithms requiring at least two qualitative variables.
    #'   - `"mixed"` for algorithms supporting both types.
    #'
    #' @return
    #' Invisibly returns `TRUE` if the dataset meets the required conditions.
    #' Throws an error otherwise.
    #'
    #' @examples
    #' \dontrun{
    #' my_algo <- MyAlgo$new()
    #' my_algo$load_and_check_data(iris)
    #' my_algo$validate_algorithm_requirements("quant")
    #' }
    #'
    #' @noRd
    validate_algorithm_requirements = function(type) {

      has_quanti <- length(private$.quanti_indices) > 1
      has_quali <- length(private$.quali_indices) > 1


      if (!type %in% c("quant", "qual", "mixed")) {
        stop("Invalid requirement type.")
      }

      if (type == "quant" && length(private$.quanti_indices) < 2) {
        stop("Not enough quantitatives columns. Consider changing vartype.")
      }

      if (type == "qual" && length(private$.quali_indices) < 2) {
        stop("Not enough qualitative columns. Consider changer vartype.")
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

    # Access or set the data attributes (original dataframe)
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

    # Access or set the number of clusters generated.
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

    # Access or modify cluster labels.
    labels = function(value) {
      if (missing(value)) {
        return(private$.labels)
      }

      # Cheking length consistency with data.
      if (!is.null(private$.data) && length(value) != ncol(private$.data)) {
        stop("'labels' length must match number of variables.")
      }

      private$.labels <- value
    },


    # -------------------------------------------------------------------------
    # fitted getter/setter
    # -------------------------------------------------------------------------

    # Access or edit fitted attribute.
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


    # access to the indices of quantitative indices.
    quanti_indices = function() private$.quanti_indices,

    # access to the indices of qualitative indices
    quali_indices  = function() private$.quali_indices
  )
)



