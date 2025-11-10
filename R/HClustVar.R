#' Hierarchical Variable Clustering
#'
#' @description
#' Implements a variable clustering algorithm using Agglomerative Hierarchical
#' Clustering (AHC). Supports quantitative, qualitative, or mixed data.
#' Quantitative variables can be clustered using correlation-based distances,
#' and qualitative variables using Cramer's V.
#'
#' @details
#' The `HClustVar` class inherits from `ClusteringBase` and provides:
#' - Automatic detection of variable types (quantitative, qualitative, mixed)
#' - Calculation of a distance matrix suitable for hierarchical clustering
#' - Internal preprocessing of quantitative variables (quantile discretization)
#' - Cramer's V calculation for categorical associations
#' - Storage of the resulting hierarchical clustering tree
#'
#' @section Public Methods:
#' \describe{
#'   \item{\code{initialize(vartype = "auto", dist.metric = NULL)}}{
#'     Creates a new instance of HClustVar.
#'     \itemize{
#'       \item \code{vartype}: Type of variables ("auto", "quant", "qual", "mixed")
#'       \item \code{dist.metric}: Distance metric ("rsquare", "r") for vartype = "quant"
#'     }
#'   }
#'   \item{\code{fit(data)}}{
#'     Fits the clustering model on the provided data.
#'     \itemize{
#'       \item \code{data}: Data frame containing the variables to cluster
#'     }
#'     This method:
#'     \itemize{
#'       \item Validates and loads the data
#'       \item Automatically determines variable type if vartype = "auto"
#'       \item Calculates the appropriate distance matrix
#'       \item Performs hierarchical clustering using Ward's method
#'     }
#'   }
#'   }
#'   \item{\code{plot_dendrogram(k = NULL, ...)}}{
#'     Displays the dendrogram of the hierarchical clustering.
#'     \itemize{
#'       \item \code{k}: Number of clusters to highlight (optional)
#'       \item \code{...}: Additional arguments passed to plot()
#'     }
#'   }
#' }
#'
#' @section Active Bindings (getters):
#' \describe{
#'   \item{\code{dist.metric}}{
#'     Returns the distance metric used ("rsquare" or "r" for quantitative variables)
#'   }
#'   \item{\code{vartype}}{
#'     Returns the type of variables processed ("quant", "qual", or "mixed")
#'   }
#'   \item{\code{dist.matrix}}{
#'     Returns the calculated distance matrix (object of class "dist")
#'   }
#' }
#'
#' @section Behavior by Variable Type:
#' \subsection{Quantitative variables (vartype = "quant")}{
#'   \itemize{
#'     \item Calculates the correlation matrix between variables
#'     \item Applies the chosen metric (r or rsquare)
#'     \item Creates a dissimilarity matrix: \code{dist = sqrt(1 - cor)}
#'   }
#' }
#' \subsection{Qualitative variables (vartype = "qual")}{
#'   \itemize{
#'     \item Calculates Cramer's V between all pairs of variables
#'     \item Creates a dissimilarity matrix: \code{dist = 1 - cramer_v}
#'   }
#' }
#' \subsection{Mixed variables (vartype = "mixed")}{
#'   \itemize{
#'     \item Discretizes quantitative variables into 4 quantiles
#'     \item Treats all variables as qualitative
#'     \item Calculates Cramer's V on the entire dataset
#'     \item Creates a dissimilarity matrix: \code{dist = 1 - cramer_v}
#'   }
#' }
#' \subsection{Automatic detection (vartype = "auto")}{
#'   \itemize{
#'     \item Selects "quant" if all variables are quantitative
#'     \item Selects "qual" if all variables are qualitative
#'     \item Selects "mixed" if both types are present
#'   }
#' }
#'
#' @section Clustering Algorithm:
#' Ward's method (ward.D) is used for hierarchical clustering.
#' This method minimizes within-cluster variance and produces
#' relatively homogeneous clusters.
#'
#' @section Error and Warning Handling:
#' \itemize{
#'   \item Error if \code{dist.metric} is invalid for vartype = "quant"
#'   \item Error if \code{vartype} is invalid
#'   \item Warning if \code{dist.metric} is specified for vartype "qual" or "mixed"
#'   \item Warning if automatic detection changes the metric
#' }
#'
#' @examples
#' \dontrun{
#' # Example with quantitative variables
#' data_quant <- data.frame(
#'   var1 = rnorm(100),
#'   var2 = rnorm(100),
#'   var3 = rnorm(100)
#' )
#'
#' hc_quant <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
#' hc_quant$fit(data_quant)
#' hc_quant$plot_dendrogram()
#'
#' # Example with qualitative variables
#' data_qual <- data.frame(
#'   cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
#'   cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
#'   cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
#' )
#'
#' hc_qual <- HClustVar$new(vartype = "qual")
#' hc_qual$fit(data_qual)
#' hc_qual$plot_dendrogram()
#'
#' # Example with automatic detection
#' data_mixed <- data.frame(
#'   num1 = rnorm(100),
#'   num2 = rnorm(100),
#'   cat1 = factor(sample(letters[1:3], 100, replace = TRUE))
#' )
#'
#' hc_auto <- HClustVar$new(vartype = "auto")
#' hc_auto$fit(data_mixed)
#' hc_auto$plot_dendrogram()
#'
#' # Access properties
#' print(hc_auto$vartype)      # "mixed"
#' print(hc_auto$dist.matrix)  # Distance matrix
#' }
#'
#' @section Dependencies:
#' This class requires:
#' \itemize{
#'   \item The R6 package
#'   \item The parent class ClusteringBase
#'   \item Base R functions: cor(), hclust(), chisq.test(), table()
#' }
#'
#' @note
#' This class is not exported; it is intended for internal use
#' within the package.
#'
#' @seealso
#' \code{\link{ClusteringBase}} for the parent class
#' \code{\link[stats]{hclust}} for hierarchical clustering
#' \code{\link[stats]{cor}} for correlation calculations
#'
#' @references
#' \itemize{
#'   \item Abdallah, H. ; Saporta, G. Revue de Statistique Appliquée, Tome 46 (1998) no. 4, pp. 5-26
#'   \item Rakotomalala R. Tutoriel Tanagra, « Classification de variables qualitatives», décembre 2013
#' }
#'
#' @family clustering classes
#' @keywords internal
#'
#' @noRd
HClustVar <- R6::R6Class(
  "HClustVar",
  inherit = ClusteringBase,

  # ==========================================================================
  # PRIVATE FIELDS
  # ==========================================================================
  private = list(

    # ==== Parameters

    # Metric to calculate distances between variables.
    .dist.metric = NULL,

    # Type of variable to make clustering on.
    .vartype = NULL,
    .auto_var_selection = FALSE,

    # CAH method chosen.
    .cah.method = NULL,

    # ==== Distances calculations

    # Contains the distances matrix.
    .dist.matrix = NULL,

    # Keep the HAC object.
    .tree = NULL,

    # centroids (computed after cut_tree function, latent components of each cluster as dataframe stored as a list of vectors).
    .centroids = NULL,

    # a vector of eigen value of the latent component of each cluster.
    .clusters.eigen = NULL,


    # ==========================================================================
    # PRIVATE METHODS
    # ==========================================================================

    # -----------------------------------------------------------------------
    # Check input method
    # -----------------------------------------------------------------------

    # Validate input parameters for class instantiation
    # Checks that 'vartype' and 'dist.metric' are valid.
    # - vartype: 'quant', 'qual', 'mixed', or 'auto'
    # - dist.metric: 'r' or 'rsquare' (only if vartype = 'quant')
    # - cah.method: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'.
    check_input = function(vartype, dist.metric, cah.method) {

      # Checking metric selection
      if (vartype == "quant" && !(dist.metric %in% c("rsquare", "r"))) {
        stop(paste0("Parameter 'dist.metric' has invalid value for quantitative variables. Choose: 'rsquare', 'r'. Got: ", dist.metric))

        # Checking vartype selection
      } else if (!(vartype %in% c("quant", "qual", "mixed", "auto"))) {
        stop(paste0("Parameter 'vartype' has invalid value. Choose: 'quant', 'qual', 'mixed', 'auto'. Got: ", vartype))

        # Warn if dist.metric is specified but ignored
      } else if (vartype %in% c("qual", "mixed") && !is.null(dist.metric)) {
        warning("'dist.metric' parameter will be ignored as it only affects 'quant' vartype")
      }

      # Test CAH method selection.
      if (!cah.method %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')) {
        stop("Unknow CAH method selected. Please select: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'")
      }
    },

    # -----------------------------------------------------------------------
    # Type auto-detection method
    # -----------------------------------------------------------------------
    # TODO: documentation
    auto_detect_vartype = function() {

      n_quanti <- length(private$.quanti_indices)
      n_quali <- length(private$.quali_indices)

      if (n_quali == 0) {

        # Only quantitative
        if (is.null(private$.dist.metric)) {
          private$.dist.metric <- "rsquare"

        } else if (!private$.dist.metric %in% c('r', 'rsquare')) {
          warning("Invalid metric for quantitative data. Using 'rsquare'.")
          private$.dist.metric <- "rsquare"
        }
        return("quant")

      } else if (n_quanti == 0) {

        # Only qualitative
        if (!is.null(private$.dist.metric)) {
          warning("dist.metric ignored for qualitative data")
        }
        return("qual")

      } else {
        # Mixte
        return("mixed")
      }
    },


    # -----------------------------------------------------------------------
    # distance_matrix computing method
    # -----------------------------------------------------------------------
    # TODO: documentation
    compute_distance_matrix = function() {

      if (private$.vartype == "quant") {
        cor_matrix <- cor(self$get_quanti_data())

        if (private$.dist.metric == "rsquare") {
          cor_matrix <- cor_matrix^2
        }

        return(as.dist(sqrt(1 - cor_matrix)))

      } else if (private$.vartype %in% c("qual", "mixed")) {

        df_quali <- if (private$.vartype == "qual") {
          self$get_quali_data()
        } else {
          private$quantile_discretisation(self$data, self$quanti_indices, 4)
        }

        vmatrix <- private$cramer_matrix(df_quali)
        return(as.dist(1 - vmatrix))

      } else {
        stop("Invalid vartype in object state")
      }
    },


    # -----------------------------------------------------------------------
    # Discretisation method
    # -----------------------------------------------------------------------

    # Transform quantitative columns into discrete variables (quantile-based). Leave unchanged one hot encoded columns.
    # - df: data.frame containing the data
    # - quanti_index: integer vector of column indices for quantitative variables
    # - n_groups: number of quantile groups to split each variable into
    # Returns a data.frame with quantitative columns discretized into factor levels.
    quantile_discretisation = function(df, quanti_index, n_groups) {

      # If no quantitative variable specified of quanti_index are not numeric, change nothing.
      if (length(quanti_index) == 0) {
        return(df)
      }

      # Create a copy of the original dataframe.
      df_copy <- as.data.frame(df)

      for (i in quanti_index) {

        # Check if the numerical value is one hot encoded.
        if (length(unique(df[[i]])) <= 2 && all(df[[i]] %in% c(0, 1))) {
          df_copy[[i]] <- as.factor(df[[i]])
          next
        }

        if (!is.numeric(df[[i]])) {
          next
        }

        # Get the quantiles.
        # Use unique to avoid small dataset creating duplicate quantile number of groups too high.
        quantiles <- unique(quantile(df[[i]], probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE))

        # transform the quantitative columns into qualitatives.
        df_copy[[i]] <- cut(
            df[[i]],

            # Set the breaks on the desired quantile of the column.
            breaks = quantiles,
            include.lowest = TRUE,
            labels = paste0("Q", 1:(length(quantiles) - 1))  # Create the new label.
          )

        # Ensure the factor conversion.
        df_copy[[i]] <- as.factor(df_copy[[i]])
      }

      return(df_copy)
    },


    # -----------------------------------------------------------------------
    # Cramer's V matrix calculation.
    # -----------------------------------------------------------------------

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

      # Check if we have enough data.
      if (length(x) == 0 || length(y) == 0) {
        warning("Empty vectors provided to cramer_v, returning 0")
        return(0)
      }

      # Create the contingency table.
      contingency <- table(x, y)

      # Check table validity
      if (nrow(contingency) < 2 || ncol(contingency) < 2) {
        warning("Contingency table too small for chi-square test, returning 0")
        return(0)
      }

      # Catch any error if computing CHISquare test fail.
      chi2_result <- tryCatch(
        suppressWarnings(chisq.test(contingency, correct = FALSE)),
        error = function(e) {
          warning(sprintf("Chi-square test failed: %s. Returning 0.", e$message))
          return(list(statistic = 0))
        }
      )

      # Get the chi2 statistic.
      chi2 <- chi2_result$statistic

      # Get the total count.
      n <- sum(contingency)

      # Get the min dimension for the dof.
      min_dim <- min(nrow(contingency), ncol(contingency)) - 1

      # Avoid dividing by 0
      if (min_dim == 0 || n == 0) {
        return(0)
      }

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

      # Initialize cramer's matrix (diag to 1)
      cramer_mat <- diag(n_vars)
      dimnames(cramer_mat) <- list(var_names, var_names)

      # Get all the combination indices (to fill up the upper triangle)
      indices <- combn(n_vars, 2)

      # Iteration over all the indices
      for (k in seq_len(ncol(indices))) {
        i <- indices[1, k]
        j <- indices[2, k]

        # Compute the cramer v
        v <- private$cramer_v(df[[i]], df[[j]])

        # update both side of the symmetric matrix.
        cramer_mat[i, j] <- v
        cramer_mat[j, i] <- v
      }
      return(cramer_mat)
    },

    # -----------------------------------------------------------------------
    # Compute centroids method.
    # -----------------------------------------------------------------------


    # TODO: Documentation
    # Compute centroids when cutting tree
    compute_centroids = function() {

      if (is.null(self$labels)) {
        stop("Cannot compute centroids: tree has not been cut yet")
      }

      # Initialization
      n_obs <- nrow(self$data)
      n_clust <- self$n_clusters

      centroids <- matrix(
        nrow = n_obs,
        ncol = n_clust,
        dimnames = list(NULL, paste0("C", 1:n_clust))
      )

      eigens <- numeric(n_clust)
      names(eigens) <- paste0("C", 1:n_clust)

      # Process each cluster
      for (i in 1:n_clust) {

        # Identify cluster variables
        cluster_vars <- names(self$labels)[self$labels == i]

        if (length(cluster_vars) == 0) {
          stop(sprintf("Cluster %d is empty", i))
        }

        # Extract cluster data
        sub_data <- self$data[, cluster_vars, drop = FALSE]

        # ═══════════════════════════════════════════════════════════════
        # Factorial analysis
        # ═══════════════════════════════════════════════════════════════
        res <- tryCatch({

          if (all(sapply(sub_data, is.numeric))) {
            FactoMineR::PCA(
              sub_data,
              ncp = 1,
              graph = FALSE,
              scale.unit = TRUE
            )

          } else if (all(sapply(sub_data, is.factor))) {
            FactoMineR::MCA(
              sub_data,
              ncp = 1,
              graph = FALSE
            )

          } else {
            FactoMineR::FAMD(
              sub_data,
              ncp = 1,
              graph = FALSE
            )
          }

        }, error = function(e) {
          stop(sprintf(
            "FactoMineR analysis failed for cluster %d: %s",
            i, e$message
          ))
        })

        # ═══════════════════════════════════════════════════════════════
        # Extract results
        # ═══════════════════════════════════════════════════════════════

        # Eigenvalue
        eigens[i] <- res$eig[1, 1]

        # First component (centroid)
        if (is.matrix(res$ind$coord)) {
          first_component <- res$ind$coord[, 1]
        } else if (is.data.frame(res$ind$coord)) {
          first_component <- res$ind$coord[[1]]
        } else {
          first_component <- as.numeric(res$ind$coord)
        }

        centroids[, i] <- first_component
      }

      # ═══════════════════════════════════════════════════════════════
      # Store results in private fields
      # ═══════════════════════════════════════════════════════════════

      private$.centroids <- centroids
      private$.clusters.eigen <- eigens

      # Do not return any value.
      invisible(NULL)
    },


    # TODO: documentation
    # Compute the correlation ratio between a qualitative variable and a quantitative variable.
    correlation_ratio = function(quali, quanti) {

      quali  <- as.factor(quali)
      quanti <- as.numeric(quanti)


      # compute total sum of squares.
      sst <- sum(
        (quanti - mean(quanti))**2
      )

      # compute sum of squares between.
      ssb <- sum(
        tapply(quanti, quali, function(v) {length(v) * (mean(v) - mean(quanti))^2})
      )

      # calculate correlation factor.
      eta2 <- ssb / sst

      return(eta2)
    }


  ),

  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list(

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------

    #' @title Constructor for HClustVar
    #'
    #' @description
    #' Initializes a new instance of the `HClustVar` class.
    #' Sets up the variable type, distance metric, and hierarchical clustering method.
    #'
    #' @param vartype Character, optional. Type of variables to cluster.
    #'   Accepted values: 'quant', 'qual', 'mixed', 'auto'. Default: 'auto'.
    #' @param dist.metric Character, optional. Distance metric for quantitative variables.
    #'   Accepted values: 'r', 'rsquare'. Ignored if vartype is not 'quant'.
    #' @param cah.method Character, optional. Method for hierarchical clustering.
    #'   Passed to `hclust()`.
    #'   Accepted values: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'
    #'   Default: 'ward.D'.
    #'
    #' @details
    #' The constructor validates input parameters using `private$check_input` and
    #' assigns them to private fields.
    #'
    #' @return A new instance of `HClustVar`.
    #' @noRd
    initialize = function(vartype = "auto", dist.metric = NULL, cah.method = "ward.D") {

      # update the auto_var_selection to allow refit with auto method.
      if (vartype == "auto") {
        private$.auto_var_selection <- TRUE
      }


      # Select default parameter for the metric if NULL.
      if (vartype == "quant" && is.null(dist.metric)) {
        dist.metric <- "rsquare"
      }

      private$check_input(vartype, dist.metric, cah.method)

      private$.dist.metric <- dist.metric
      private$.vartype <- vartype
      private$.cah.method <- cah.method
    },

    # -----------------------------------------------------------------------
    # Fit method
    # -----------------------------------------------------------------------

    #' Fit hierarchical clustering on variables
    #'
    #' @description
    #' Computes the hierarchical clustering based on the variable type (quantitative,
    #' qualitative, or mixed) and the specified distance metric.
    #'
    #' @param data data.frame or matrix containing the variables to cluster. Be careful to set your qualitatives or dummified data as factors.
    #'
    #' @details
    #' - If `vartype = "auto"`, the method automatically selects 'quant', 'qual', or 'mixed'
    #'   based on the input data.
    #' - For quantitative variables, computes correlation or squared correlation as dissimilarity.
    #' - For qualitative or mixed variables, transforms quantitative columns into quantiles (discretization)
    #'   and computes Cramer's V dissimilarity matrix.
    #' - Finally, applies hierarchical clustering (`hclust`) using the method defined in `.cah.method`.
    #'
    #' @return None. The clustering tree is stored internally in `private$.tree`.
    #' @noRd
    fit = function(data) {

      self$load_and_check_data(data)

      # vartype autodetection
      if (private$.auto_var_selection == TRUE) {
        private$.vartype <- private$auto_detect_vartype()
      }

      # Validation
      self$validate_algorithm_requirements(private$.vartype)

      # Metric distance calculation
      private$.dist.matrix <- private$compute_distance_matrix()

      # CAH
      private$.tree <- hclust(private$.dist.matrix, method = private$.cah.method)
      self$fitted <- TRUE
    },


    # -----------------------------------------------------------------------
    # cut_tree method
    # -----------------------------------------------------------------------

    #' Cut the hierarchical clustering tree into clusters
    #'
    #' @description
    #' Cuts the hierarchical clustering tree of a fitted model into clusters and compute their centroids.
    #' This method updates the `labels` and `n_clusters`attributes of the object with cluster assignments.
    #'
    #' @details
    #' The method uses `cutree()` internally to generate cluster labels for each observation.
    #' You can specify either the number of clusters (`k`) or a cut height (`h`):
    #' - If `k` is provided, the tree is cut to produce exactly `k` clusters.
    #' - If `h` is provided, the tree is cut at the specified height, and the number of clusters
    #'   is determined automatically.
    #' - If both `k` and `h` are provided, `cutree()` will prioritize `k`.
    #'
    #' The model must be fitted beforehand (`self$fitted` must be TRUE), otherwise an error
    #' is raised.
    #'
    #' @param k Integer, optional. Number of clusters to generate.
    #' @param h Numeric, optional. Height at which to cut the tree.
    #'
    #' @return A vector of cluster labels corresponding to each observation. The labels are also
    #'         stored in the object's `labels` attribute.
    #'
    #' @examples
    #' \dontrun{
    #' # Using an instance of the class:
    #' obj$cut_tree(k = 3)       # cut into 3 clusters
    #' obj$cut_tree(h = 150)     # cut at height 150
    #' }
    #'
    #' @noRd
    cut_tree = function(k = NULL, h = NULL) {

      # Check if the model is fitted.
      if (!self$fitted) {
        stop("Your model should be fitted on data first.")
      }

      # Update the labels attribute.
      self$labels <- cutree(self$tree, k, h)

      # update the number of clusters
      self$n_clusters <- length(unique(self$labels))

      # Compute the centroids of the clusters.
      private$compute_centroids()

      return(self$labels)

    },

    # -----------------------------------------------------------------------
    # plot dendrogram method
    # -----------------------------------------------------------------------

    #' Plot dendrogram of the hierarchical clustering
    #'
    #' @description
    #' Displays the dendrogram of the computed hierarchical clustering tree,
    #' with optional colored rectangles highlighting the cluster groups.
    #'
    #' @details
    #' The dendrogram visualizes the hierarchical relationships between variables
    #' obtained from the hierarchical clustering algorithm.
    #' If a number of clusters `k` is provided, the function colors and outlines
    #' the corresponding groups in the dendrogram using distinct colors.
    #'
    #' @param k (integer) Optional. Number of clusters to highlight in the dendrogram.
    #'        If `k = 0` (default), the dendrogram is drawn without colored groups.
    #'
    #' @return None. Generates a plot.
    #'
    #' @examples
    #' # Example usage (assuming the model has been fitted):
    #' # obj$plot_dendrogram(k = 3)
    #'
    #' @noRd
    plot_dendrogram = function(k = 0) {
      plot(private$.tree, main = "Dendrogramme - CAH", xlab = "Variables")

      # If k is specified, draw the rectangles.
      if (k > 0 && is.numeric(k)) {
        rect.hclust(private$.tree, k = as.integer(k), border = rainbow(k))
      }
    },

    # -----------------------------------------------------------------------
    # Predict method
    # -----------------------------------------------------------------------

    #' Predict cluster membership for new variables
    #'
    #' @description
    #' Assigns new variables to existing clusters based on their similarity
    #' with the variables used to build the clustering model.
    #'
    #' @details
    #' The prediction method varies depending on the clustering method used:
    #' - For Ward methods ('ward.D', 'ward.D2', 'centroid'): Uses correlation
    #'   with cluster centroids (principal components)
    #' - For other methods ('single', 'complete', 'median', 'average', 'mcquitty'):
    #'   Uses correlation/Cramer's V with individual variables in each cluster
    #'
    #' The similarity metric depends on variable types:
    #' - Quantitative vs Quantitative: Pearson correlation (r or r²)
    #' - Qualitative vs Qualitative: Cramer's V
    #' - Quantitative vs Qualitative: Correlation ratio (η²)
    #'
    #' @param new_data A data.frame or matrix containing new variables to classify.
    #'   Must have the same number of rows as the training data.
    #'
    #' @return A named numeric vector where:
    #'   - Names correspond to column names in new_data
    #'   - Values are cluster assignments (integers from 1 to n_clusters)
    #'
    #' @examples
    #' \dontrun{
    #' # Fit model
    #' hc <- HClustVar$new(vartype = "quant")
    #' hc$fit(training_data)
    #' hc$cut_tree(k = 3)
    #'
    #' # Predict new variables
    #' new_vars <- data.frame(new_v1 = rnorm(100), new_v2 = rnorm(100))
    #' predictions <- hc$predict(new_vars)
    #' print(predictions)
    #' # new_v1 new_v2
    #' #      1      2
    #' }
    #'
    #' @noRd
    predict = function(new_data) {

      # ==========================================================================
      # 1. Prerequisites validation
      # ==========================================================================

      # Check if model is fitted
      if (!self$fitted) {
        stop("Model must be fitted before prediction. Use fit() first.")
      }

      # Check if labels exist (tree has been cut)
      if (is.null(self$labels)) {
        stop("Tree must be cut before prediction. Use cut_tree() first.")
      }

      # Validate new_data type
      if (!is.data.frame(new_data) && !is.matrix(new_data)) {
        stop("new_data must be a data.frame or matrix")
      }

      # Convert matrix to data.frame if necessary
      if (is.matrix(new_data)) {
        new_data <- as.data.frame(new_data)
      }

      # Check row compatibility
      if (nrow(new_data) != nrow(self$data)) {
        stop(sprintf(
          "new_data must have the same number of rows as training data (%d rows expected, got %d)",
          nrow(self$data),
          nrow(new_data)
        ))
      }

      # Check for empty new_data
      if (ncol(new_data) == 0) {
        stop("new_data must contain at least one variable")
      }

      # ==========================================================================
      # 2. INITIALIZATION
      # ==========================================================================

      result <- integer(ncol(new_data))
      names(result) <- colnames(new_data)


      # If there's only one cluster, all variables belong to it, then return.
      if (self$n_clusters == 1) {
        result[] <- 1
        return(result)
      }

      # ==========================================================================
      # 3. PREDICTION DEPENDING OF CAH METHOD
      # ==========================================================================

      # -----------
      # 3A. CENTROIDS BASED METHODS (Ward, Centroid)
      # -----------
      if (private$.cah.method %in% c("ward.D", "ward.D2", "centroid")) {

        # Validate that centroids exist
        if (is.null(self$centroids)) {
          stop("Centroids not computed. This should not happen. Object may be broken.")
        }

        # Process each new variable
        for (i in seq_len(ncol(new_data))) {

          var_name <- colnames(new_data)[i]

          # --- Case 1: New variable is quantitative ---
          if (is.numeric(new_data[[i]])) {

            # Calculate correlation with each centroid
            correlations <- apply(self$centroids, 2, function(centroid) {
              cor_val <- cor(centroid, new_data[[i]], use = "complete.obs")

              # Handle NA values
              if (is.na(cor_val)) {
                warning(sprintf(
                  "Correlation is NA for variable '%s'. Using 0 as fallback.",
                  var_name
                ))
                return(0)
              }

              return(cor_val)
            })

            # Apply metric transformation if needed
            if (!is.null(self$dist.metric) && self$dist.metric == "rsquare") {
              correlations <- correlations^2
            }

            # Assign to cluster with maximum correlation
            result[i] <- which.max(correlations)

            # --- Case 2: New variable is qualitative ---
          } else {

            # Calculate correlation ratio with each centroid
            cor_ratios <- apply(self$centroids, 2, function(centroid) {
              private$correlation_ratio(new_data[[i]], centroid)
            })

            # Assign to cluster with maximum correlation ratio
            result[i] <- which.max(cor_ratios)
          }
        }

        # ---------------
        # 3B. OTHER METHODS (Single, Complete, Median, Average, McQuitty)
        # ---------------
      } else if (private$.cah.method %in% c("single", "complete", "median", "average", "mcquitty")) {

        # Process each new variable
        for (i in seq_len(ncol(new_data))) {

          var_name <- colnames(new_data)[i]
          similarity <- numeric(ncol(self$data))
          names(similarity) <- colnames(self$data)

          # --------------
          # Calculate similarity between new variable and all training variables
          # --------------

          # --- Scenario 1: Training data is purely QUANTITATIVE ---
          if (self$vartype == "quant") {

            if (is.numeric(new_data[[i]])) {
              # Quantitative vs Quantitative: Use correlation
              similarity <- apply(self$data, 2, function(x) {
                cor(new_data[[i]], x, use = "complete.obs")
              })

              # Apply metric transformation
              if (self$dist.metric == "rsquare") {
                similarity <- similarity^2
              }

            } else {
              # Qualitative vs Quantitative: Use correlation ratio
              similarity <- apply(self$data, 2, function(x) {
                private$correlation_ratio(new_data[[i]], x)
              })
            }

            # --- Scenario 2: Training data is purely QUALITATIVE ---
          } else if (self$vartype == "qual") {

            if (is.numeric(new_data[[i]])) {
              # Quantitative vs Qualitative: Use correlation ratio (reversed)
              similarity <- apply(self$data, 2, function(x) {
                private$correlation_ratio(x, new_data[[i]])
              })

            } else {
              # Qualitative vs Qualitative: Use Cramer's V
              similarity <- apply(self$data, 2, function(x) {
                private$cramer_v(x, new_data[[i]])
              })
            }

            # --- Scenario 3: Training data is MIXED ---
          } else if (self$vartype == "mixed") {

            # Discretize new variable if quantitative
            if (is.numeric(new_data[[i]])) {
              new_var_discretized <- private$quantile_discretisation(
                df = data.frame(temp_var = new_data[[i]]),
                quanti_index = 1,
                n_groups = 4
              )[[1]]
            } else {
              new_var_discretized <- new_data[[i]]
            }

            # Calculate Cramer's V with all training variables (after discretization)
            similarity <- sapply(seq_len(ncol(self$data)), function(j) {

              train_var <- self$data[[j]]

              # Discretize training variable if quantitative
              if (is.numeric(train_var)) {
                train_var_discretized <- private$quantile_discretisation(
                  df = data.frame(temp_var = train_var),
                  quanti_index = 1,
                  n_groups = 4
                )[[1]]
              } else {
                train_var_discretized <- train_var
              }

              # Calculate Cramer's V
              private$cramer_v(train_var_discretized, new_var_discretized)
            })

            names(similarity) <- colnames(self$data)

          } else {
            stop(sprintf("Unknown vartype: %s", self$vartype))
          }

          # ----------
          # Aggregate similarities by cluster according to CAH method
          # ------------

          cluster_proximity <- numeric(self$n_clusters)

          for (cluster_id in seq_len(self$n_clusters)) {

            # Get similarities for variables in this cluster
            cluster_similarities <- similarity[self$labels == cluster_id]

            # Aggregate according to method
            cluster_proximity[cluster_id] <- switch(
              private$.cah.method,
              "single"   = max(cluster_similarities),      # Closest element
              "complete" = min(cluster_similarities),      # Farthest element
              "median"   = median(cluster_similarities),   # Median similarity
              "average"  = mean(cluster_similarities),     # Mean similarity
              "mcquitty" = mean(cluster_similarities),     # Same as average (no ponderation for illustrative var).
              stop(sprintf("Unhandled CAH method: %s", private$.cah.method))
            )
          }

          # Assign to cluster with maximum proximity
          result[i] <- which.max(cluster_proximity)
        }

        # ---------
        # 3C. UNKNOWN METHOD
        # ----------
      } else {
        stop(sprintf(
          "Unknown or unsupported CAH method: %s. Supported methods are: %s",
          private$.cah.method,
          paste(c("ward.D", "ward.D2", "centroid", "single", "complete",
                  "median", "average", "mcquitty"), collapse = ", ")
        ))
      }

      return(result)
    },

    # -----------------------------------------------------------------------
    # Within inertia plot
    # -----------------------------------------------------------------------

    # TODO: Documentation
    plot_cohesion = function() {

      # Check if fitted.
      if (!self$fitted) {
        stop("Model should be fitted on data.")
      }

      # Save the original class configuration to reset after plotting.
      original_label <- self$labels
      original_n_clusters <- self$n_clusters

      original_centroids <- self$centroids
      original_clusters.eigen <- self$clusters.eigen

      # initialize the within inertia vector.
      within.inertia <- numeric()

      max_clusters <- 0

      # loop through each cluster.
      if (self$vartype == "quant"){
        max_clusters <- length(self$quanti_indices)
      } else if (self$vartype == "qual") {
        max_clusters <- length(self$quali_indices)
      } else {
        max_clusters <- ncol(self$data)
      }

      # initialize the inertia vector
      inertia <- numeric(max_clusters)

      for (k in seq_len(max_clusters)) {


        # if we reach max cluster, just update knowing each eigens values is equal 1.
        if (k == max_clusters) {
          inertia[k] <- k
        } else {

          # cut the tree with k clusters.
          self$cut_tree(k)

          # save the sum of within inertia (eigens values of each cluster).
          inertia[k] <- sum(self$clusters.eigen)
        }
      } # end for.

      # calculate the cohesion gain.
      cohesion <- (inertia - inertia[1]) / (max_clusters - inertia[1])

      # plot graph
      plot(cohesion, type = "b", xlab = "Number of clusters", ylab = "Gain in cohesion (%)")

      # Restore original configuration.
      self$labels <- original_label
      self$n_clusters <- original_n_clusters

      self$centroids <- original_centroids
      self$clusters.eigen <- original_clusters.eigen

    },


    # -----------------------------------------------------------------------
    # Summary method
    # -----------------------------------------------------------------------

    # TODO: documentation
    summary = function() {
      # TODO: Add more element to the summary function.

      # Model state
      if (!self$fitted) {
        cat("Status: NOT FITTED\n")
        cat("Use $fit(data) to train the model.\n\n")
        return(invisible(self))
      }

      if (is.null(self$n_clusters)) {
        cat("Status: Tree not cut.\n")
        cat("Use self$cut_tree().\n\n")
        return(invisible(self))
      }

      if (self$vartype == "quant") {
        variables <- colnames(self$data[self$quanti_indices])
      } else if (self$vartype == "qual") {
        variables <- colnames(self$data[self$quali_indices])
      } else {
        variables <- colnames(self$data)
      }


      # ----------------
      # Cluster summary
      # ----------------

      # Compute the cluster summary.
      clust_summary <- data.frame(
        cluster = 1:self$n_clusters,
        n_members = sapply(1:self$n_clusters, FUN = function(x) {sum(self$labels == x)}),
        var_explained = round(self$clusters.eigen, 2)
      )

      clust_summary["prop_explained"] <- round(clust_summary$var_explained / clust_summary$n_members, 2)

      total_var_explained <- sum(clust_summary$var_explained)
      percentage_var_explained <- total_var_explained / length(variables)


      # ----------------
      # Cluster Members
      # ----------------

      # Compute Cluster members
      clust_members <- data.frame(
        member = variables,
        cluster = self$labels[variables]
      )

      own_cluster_R2 <- numeric(length(variables))
      next_closest_R2 <- numeric(length(variables))


      # Iteration on all variables to calculate theior properties.
      for (i in 1:length(variables)) {

        variable <- variables[i]

        # initialize the vector of correlation with the other clusters
        cor_clusters_temp = numeric(self$n_clusters)

        # own_cluster
        own_clust <- self$labels[variable]

        # Compute the correlation² with all the other centroids
        # Correlation if quantitative.
        if (is.numeric(self$data[[variable]])) {  # Double crochet [[]]
          cor_clusters_temp <- apply(self$centroids, 2, function(centroid) {
            cor(self$data[[variable]], centroid, use = "complete.obs")^2
          })
        } else {
          # correlation ratio if qualitative.
          cor_clusters_temp <- apply(self$centroids, MARGIN = 2, FUN = function(centroid) {private$correlation_ratio(self$data[[variable]], centroid)})
        }

        print(cor_clusters_temp)

        own_cluster_R2[i] <- cor_clusters_temp[own_clust]
        next_closest_R2[i] <- max(cor_clusters_temp[-own_clust])
      } # end for

      # Update the dataframe.
      clust_members["own_cluster_R2"] <- round(own_cluster_R2, 2)
      clust_members["next_closest_R2"] <- round(next_closest_R2, 2)
      clust_members["1 - R2_ratio"] <- round((1 - clust_members["own_cluster_R2"]) / (1 - clust_members["next_closest_R2"]), 2)

      # Sort the clust_members df
      clust_members <- clust_members[order(clust_members$cluster, -clust_members$own_cluster_R2), ]

      return(list(clust_summary = clust_summary, clust_members = clust_members))
    },

    # -----------------------------------------------------------------------
    # Print method
    # -----------------------------------------------------------------------


    #' Print method for HClustVar objects
    #'
    #' @description
    #' Displays a concise overview of the HClustVar object with key information.
    #'
    #' @return Invisibly returns the object itself.
    #'
    #' @noRd
    print = function() {

      cat("\n")
      cat("══════════════════════════════════════════════════════\n")
      cat("            Hierarchical Variable Clustering          \n")
      cat("══════════════════════════════════════════════════════\n\n")

      # Model state
      if (!self$fitted) {
        cat("Status: NOT FITTED\n")
        cat("Use $fit(data) to train the model.\n\n")
        return(invisible(self))
      }

      cat("Status: ✓ FITTED\n\n")

      # Base info
      cat(sprintf("Data:          %d variables × %d observations\n",
                  ncol(self$data), nrow(self$data)))
      cat(sprintf("Variable type: %s\n", private$.vartype))
      cat(sprintf("CAH method:    %s\n", private$.cah.method))

      if (!is.null(private$.dist.metric)) {
        cat(sprintf("Distance:      %s\n", private$.dist.metric))
      }

      # Cluster info
      if (!is.null(self$n_clusters) && !is.null(self$labels)) {
        cat(sprintf("\nClusters:      %d clusters created\n", self$n_clusters))

        # Display distribution
        cluster_sizes <- table(self$labels)
        size_str <- paste(sprintf("%d", cluster_sizes), collapse = " | ")
        cat(sprintf("Distribution:  %s variables\n", size_str))

      } else {
        cat("\nClusters:      Tree not cut (use $cut_tree())\n")
      }

      cat("\n══════════════════════════════════════════════════════\n")
      cat("Use $summary() for detailed information\n")
      cat("Use $plot_dendrogram() to visualize the tree\n\n")

      invisible(self)
    }
  ),

  # ==========================================================================
  # ACTIVE BINDINGS (getters/setters)
  # ==========================================================================

  active = list(

    # Returns the distance metric used for quantitative variables ('r' or 'rsquare')
    dist.metric = function() {return(private$.dist.metric)},

    # Returns the type of variables considered for clustering ('quant', 'qual', 'mixed', 'auto')
    vartype = function() {return(private$.vartype)},

    # Returns the computed dissimilarity/distance matrix
    dist.matrix = function() {return(private$.dist.matrix)},

    # return the tree computed from the clustering.
    tree = function() {return(private$.tree)},

    # Return the method used for CAH
    cah.method = function() {return(private$.cah.method)},

    # centroids setter / getter.
    centroids = function(value) {
      if (missing(value)) {
        return(private$.centroids)
      }

      # check type
      if (!(is.data.frame(value) || is.matrix(value))) {
        stop("Centroids should be a dataframe/matrix composed of 1st principal component of each clusters.")
      }

      # Check the 'value' dataframe dimensions
      if (!nrow(value) == nrow(self$data) || !ncol(value) == self$n_clusters) {
        stop("Centroids dataframe/matrix should have the same number of row than the fitted dataframe and 1 column by clusters")
      }

      private$.centroids <- value

    },

    # centroids setter / getter.
    clusters.eigen = function(value) {
      if (missing(value)) {
        return(private$.clusters.eigen)
      }

      # check type
      if (!is.numeric(value)) {
        stop("clusters.eigen vector should contains only numerical values.")
      }

      # Check the 'value' vector length.
      if (length(value) != self$n_clusters) {
        stop("The clusters.eigen vector should contains one value per cluster.")
      }

      private$.clusters.eigen <- value
    }

  )
)
