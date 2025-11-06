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
    # CAH method chosen.
    .cah.method = NULL,
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
      if (!cah.method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) {
        stop("Unknow CAH method selected. Please select: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'")
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

      # If no quantitative variable specified.
      if (length(quanti_index) == 0) {
        return(df)
      }

      # Create a copy of the original dataframe.
      df_copy <- df

      for (i in quanti_index) {

        # Check if the numerical value is one hot encoded.
        if (length(unique(df[[i]])) <= 2 && all(df[[i]] %in% c(0, 1))) {

          # Transform into factor.
          df_copy[[i]] <- as.factor(df[[i]])

          # Ignore the loop and go to the next column.
          next

        }

        # Get the quantiles.
        # Use unique to avoid small dataset creating duplicate quantile chen number of groups too high.
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

      # Create the contingency table.
      contingency <- table(x, y)

      # Get the chi2 statistic.
      chi2 <- suppressWarnings(chisq.test(contingency, correct = FALSE)$statistic)

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
    },

    # -----------------------------------------------------------------------
    # Compute centroids method.
    # -----------------------------------------------------------------------


    # TODO: Documentation
    # Compute centroids when cutting tree
    compute_centroids = function() {

      # Initialise centroids matrix.
      centroids <- matrix(nrow = nrow(self$data), ncol = 0)

      # Compute
      eigens <- c()

      # Iteration on all clusters.
      for (i in 1:self$n_clusters) {

        cluster_vars <- names(self$labels)[self$labels == i]
        sub_data <- self$data[, cluster_vars, drop = FALSE]

        if (all(sapply(sub_data, is.numeric))) {
          res <- FactoMineR::PCA(sub_data, ncp = 1, graph = FALSE, scale.unit = TRUE)

        } else if (all(sapply(sub_data, is.factor))) {
          res <- FactoMineR::MCA(sub_data, ncp = 1, graph = FALSE)

        } else {
          res <- FactoMineR::FAMD(sub_data, ncp = 1, graph = FALSE)
        }

        # Append the result to centroids list and eigens vector.
        eigens <- c(eigens, res$eig[1, 1])
        centroids <- cbind(centroids, res$ind$coord)

        } # endfor

      # Update the colname for centroids for clarity
      colnames(centroids) <- paste0("C", 1:self$n_clusters)

      # Set the private attributes.
      self$centroids <- centroids
      self$clusters.eigen <- eigens
    },


    # TODO: documentation
    # Compute the correlation ratio between a qualitative variable and a quantitative variable.
    correlation_ratio = function(quali, quanti) {

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

      # Validate data.
      self$load_and_check_data(data)

      # if "auto" selected, choose the best vartypes.
      if (private$.vartype == "auto") {

        # Select quantitative if only quantitatives columns.
        if (length(private$.quali_indices) == 0) {
          private$.vartype = "quant"

          if (is.null(private$.dist.metric)) {
            private$.dist.metric <- "rsquare"
          } else if (!private$.dist.metric %in% c('r', 'rsquare')) {
            warning("Chosen metric does not match quantitatives data. Changing for 'rsquare'.")
            dist.metric <- "rsquare"
          }

          # Select qualitative variable.
        } else if (length(private$.quanti_indices) == 0) {
          private$.vartype <-  "qual"

          # Raise a warning if dist.metric specified.
          if (!is.null(private$.dist.metric)) {
            warning("Dataframe has only qualitative values, param. dist.metric is ignored.")
          }


        } else {
          private$.vartype <-  "mixed"
        }
      }

      # Check if selected the data contains at least 2 columns of the selected vartype.
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
      private$.tree <- hclust(private$.dist.matrix, method = private$.cah.method)

      # Update the fitted variable.
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
        stop("Your model is not fitted with any data!")
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

    # TODO: documlentation
    predict = function(new_data) {

      # --- 1. Check prerequisites ---
      # Check if model is fitted
      if (!self$fitted) {
        stop("Model must be fitted before prediction")
      }

      # Check if labels exist
      if (is.null(self$labels)) {
        stop("Tree must be cut before prediction. Use cut_tree() first.")
      }

      # Validate new_data
      if (!is.data.frame(new_data) && !is.matrix(new_data)) {
        stop("new_data must be a data.frame or matrix")
      }

      # Check row compatibility
      if (nrow(new_data) != nrow(self$data)) {
        stop("new_data must have the same number of rows as training data")
      }

      # For each variable of the new data, compute the correlation for each clusters and assign them to the closest.

      # --- 2. Initialize result vector ---
      result <- c()


      # --- 3. Iterate over each variable (column) in the new dataset ---
      # Iteration on each new variable.
      for (i in 1:ncol(new_data)) {

        # --- 4. Quantitative illustrative variable --
        if (is.numeric(new_data[, i])) {

          # Check the dist.metric
          if (!is.null(self$dist.metric) && self$dist.metric == "r") {
            correlations <- apply(self$centroids, 2, FUN = function(x) {cor(x, new_data[, i])})

            # add the maximum index to the result vector.
            result <- c(result, which.max(correlations))

          } else {
            correlations <- apply(self$centroids, 2, FUN = function(x) {cor(x, new_data[, i])**2})


            # add the maximum index to the result vector.
            result <- c(result, which.max(correlations))
          } # end if


          # --- 5. Qualitative illustrative variable ---
        } else {

          correlations <- apply(self$centroids, 2, FUN = function(x) {private$correlation_ratio(new_data[, i], x)})

          result <- c(result, which.max(correlations))

        }
      } # end for

      # Rename the vector of label.
      names(result) <- names(new_data)

      # Return the vector of label.
      return(result)

    } # end predict
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
