
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
    .dist.metric = NULL,

    # Type of variable to make clustering on.
    .vartype = NULL,

    # Contains the distances matrix.
    .dist.matrix = NULL,

    # Contains the method to calculate the distance matrix.
    .dist.method = NULL,

    # Keep the HAC object.
    .tree = NULL,


    # ==========================================================================
    # PRIVATE MEHTODS
    # ==========================================================================

    #' Validate input parameters for class instanciation.
    #' @param vartype The type of variable to make the clustering on. Accepted values: ['quant', 'qual', 'mixed'].
    #' @param dist.metric The metric to use on the algorithm. Only available when vartype = 'quant'. Accepted values:
    check_input = function(vartype, dist.metric) {

      # Checking metric selection
      if (vartype == "quant" && !(dist.metric %in% c("rsquare", "r"))) {
        stop(paste0("Parameter 'dist.metric' has invalid value. Choose: 'rsquare', 'r'. Got: ", dist.metric))

        # Checking vartype selection
      } else if (!(vartype %in% c("quant", "qual", "mixed"))) {
        stop(paste0("Parameter 'vartype' has invalid value. Choose: 'quant', 'qual', 'mixed'. Got: ", vartype))

        # Warn if dist.metric is specified but ignored
      } else if (vartype %in% c("qual", "mixed") && !is.null(dist.metric)) {
        warning("'dist.metric' parameter will be ignored as it only affects 'quant' vartype")
      }
    },

    #' Transform the quantitative columns of a dataframe into a discrete variable.
    #' Split the variable into 4 quantiles.
    quantile_discretisation = function(df, quanti_index, n_groups) {

      # If no quantitative variable specified.
      if (length(quanti_index) == 0) {
        return(df)
      }

      # Create a copy of the original dataframe.
      df_copy <- df

      for (i in quanti_index) {

        # Get the quantiles.
        # Use unique to avoid small dataset creating duplicate quantile chen number of groups too high.
        quantiles <- unique(quantile(df[[i]], probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE))

        # transform the quantitative columns into qualitatives.
        df_copy[i] <- cut(
            df[[i]],

            # Set the breaks on the desired quantile of the column.
            breaks = quantiles,
            include.lowest = TRUE,
            labels = paste0("Q", 1:(length(quantiles) - 1))  # Create the new label.
          )
      }

      return(df_copy)
    },

    #' Calculate Cramer's V between two variables
    #'
    #' Internal helper function to compute Cramer's V,
    #' a measure of association between two categorical variables.
    #'
    #' @param x A factor or vector representing the first categorical variable.
    #' @param y A factor or vector representing the second categorical variable.
    #'
    #' @return A numeric value representing Cramer's V (between 0 and 1).
    cramer_v = function(x, y) {

      # Create the contingency table.
      contingency <- table(x, y)

      # Get the chi2 statistic.
      chi2 <- chisq.test(contingency, correct = FALSE)$statistic

      # Get the total count.
      n <- sum(contingency)

      # Get the min dimension for the dof.
      min_dim <- min(nrow(contingency), ncol(contingency)) - 1

      # Calculate the cramer V.
      v <- sqrt(chi2 / (n * min_dim))
      return(as.numeric(v))
    },

    #' Calculate Cramer's V matrix
    #'
    #' Internal helper method that computes a matrix of pairwise Cramer's V
    #' coefficients between all categorical variables in a data frame.
    #'
    #' @param df A data frame containing categorical variables.
    #'
    #' @return A symmetric numeric matrix where each element represents the
    #' Cramer's V value between two variables (values range from 0 to 1).
    cramer_matrix = function(df) {
      n_vars <- ncol(df)
      var_names <- colnames(df)

      cramer_mat <- matrix(1, nrow = n_vars, ncol = n_vars,
                           dimnames = list(var_names, var_names))

      for (i in 1:(n_vars - 1)) {
        for (j in (i + 1):n_vars) {
          v <- private$cramer_v(df[[i]], df[[j]])
          cramer_mat[i, j] <- v
          cramer_mat[j, i] <- v
        }
      }
      return(cramer_mat)
    }

  ),

  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list(

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------
    initialize = function(vartype = "mixed", dist.metric = NULL) {

      # Select default parameter for the metric if NULL.
      if (vartype == "quant" && is.null(dist.metric)) {
        dist.metric <- "rsquare"
      }

      private$check_input(vartype, dist.metric)

      private$.dist.metric <- dist.metric
      private$.vartype <- vartype
    },

    # -----------------------------------------------------------------------
    # Fit method
    # -----------------------------------------------------------------------
    fit = function(data) {

      # Validate data.
      self$load_and_check_data(data)

      # Check if enough var types for the
      self$validate_algorithm_requirements(private$.vartype)

      # Create the distances matrix
      if (private$.vartype == "quant") {

        # Generate the correlation matrix based on correlation.
        cor_matrix <- cor(self$get_quanti_data())

        # Use the squared correlation matrix if requested.
        if (private$.dist.metric == "rsquare") cor_matrix <- cor_matrix^2

        # Create the dissimilarity matrix used as base distances matrix.
        private$.dist.matrix <- as.dist(sqrt(1 - cor_matrix))

      } else if (private$.vartype %in% c("qual", "mixed")) {

        # Prepare the data depending if we want the mixed clustering or not.
        if (private$.vartype == "qual") {

          # Get the qualtative data only.
          df_quali <- self$get_quali_data()

        } else if (private$.vartype == "mixed") {

          # Preparation of the data by transforming quantitatives into qualitatives.
          df_quali <- private$quantile_discretisation(data, self$quanti_indices, 4)
        }

        # Compute the cramer's V matrix.
        vmatrix <- private$cramer_matrix(df_quali)

        # Generetae the dissimilarity matrix.
        private$.dist.matrix <- as.dist(1 - vmatrix)

      # Raise an error if no corresponding type.
      } else {
        stop("Invalid class structure: an error occured or one of the attribute has been manually edited.")
      }

      # Compute CAH.
      private$.tree <- hclust(private$.dist.matrix, method = "ward.D")
    },

    cut_tree = function(k = NULL, h = NULL) {
      # TODO: A implémenter
    },

    # Dendrogramme
    plot_dendrogram = function(k = NULL, ...) {
      plot(private$.tree)
    }
  ),

  # ==========================================================================
  # ACTIVE BINDINGS (getters/setters)
  # ==========================================================================

  active = list(

    dist.metric = function() {return(private$.dist.metric)},
    vartype = function() {return(private$.vartype)},
    dist.matrix = function() {return(private$.dist.matrix)}

  )
)
