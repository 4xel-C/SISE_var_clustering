#' @title Hierarchical Variable Clustering
#'
#' @description
#' Implements hierarchical clustering of variables using Agglomerative Hierarchical
#' Clustering (AHC). Groups variables based on their similarity from correlation
#' or cramer v, supporting quantitative, qualitative, and mixed data types.
#'
#' @details
#' The HClustVar class provides comprehensive tools for variable clustering with
#' automatic type detection, multiple distance metrics, and extensive diagnostics.
#'
#' Typical workflow:
#' \enumerate{
#'   \item Initialize with `new()`
#'   \item Fit with `fit(data)`
#'   \item Select k using `plot_agg_levels()`
#'   \item Cut with `cut_tree(k)`
#'   \item Analyze with `summary()` and `plot_silhouette()`
#' }
#'
#' For quantitative variables, distance metrics "r" or "rsquare" are available.
#' For qualitative variables, Cramer's V is used. For mixed data, quantitative
#' variables are discretized into quartiles before applying Cramer's V.
#'
#' Available clustering methods: ward.D (default, recommended), ward.D2, single,
#' complete, average, mcquitty, median, centroid.
#'
#' @section Workflow and Usage:
#' **Basic Workflow:**
#' \preformatted{
#' # 1. Initialize
#' hc <- HClustVar$new(vartype = "auto", dist.metric = "rsquare")
#'
#' # 2. Fit on data
#' hc$fit(my_data)
#'
#' # 3. Visualize dendrogram
#' hc$plot_dendrogram()
#'
#' # 4. Choose optimal number of clusters
#' hc$plot_agg_levels()
#'
#' # 5. Cut tree
#' hc$cut_tree(k = 3)
#'
#' # 6. Analyze results
#' results <- hc$summary()
#' hc$plot_silhouette()
#' hc$mds_projection()
#' }
#'
#' @section Variable Types:
#' The class supports four variable type configurations:
#' \describe{
#'   \item{`"auto"`}{Automatically detects variable types in the data. Chooses "quant"
#'     if all variables are numeric, "qual" if all are factors, or "mixed" if both
#'     types are present.}
#'   \item{`"quant"`}{For purely quantitative (numeric) variables. Uses correlation-based
#'     distances (Pearson correlation or squared correlation).}
#'   \item{`"qual"`}{For purely qualitative (categorical) variables. Uses Cramer's V
#'     as association measure.}
#'   \item{`"mixed"`}{For datasets containing both quantitative and qualitative variables.
#'     Discretizes numeric variables into quartiles and applies Cramer's V to all variables.}
#' }
#'
#' @section Distance Metrics:
#' Distance metrics determine how similarity between variables is measured:
#' \describe{
#'   \item{`"r"`}{Absolute Pearson correlation coefficient. Distance = sqrt(1 - r).
#'     Suitable when you want to group variables with strong linear relationships
#'     regardless of direction.}
#'   \item{`"rsquare"`}{Squared Pearson correlation (R²). Distance = sqrt(1 - r²).
#'     Default for quantitative variables. Treats positive and negative correlations
#'     equally, focusing on strength of linear relationship.}
#' }
#'
#' @section Clustering Methods:
#' The class supports all standard hierarchical clustering linkage methods:
#' \describe{
#'   \item{`"ward.D"`}{**Recommended default.** Minimizes within-cluster variance using
#'     Ward's minimum variance method (original algorithm). Tends to create compact,
#'     spherical clusters of similar size.}
#'   \item{`"ward.D2"`}{Ward's method with squared dissimilarities. Equivalent to ward.D
#'     when using Euclidean distances. May produce slightly different results with
#'     correlation-based distances.}
#'   \item{`"single"`}{Single linkage (nearest neighbor). Distance between clusters is
#'     the minimum distance between any two points. Can produce elongated "chain-like"
#'     clusters.}
#'   \item{`"complete"`}{Complete linkage (farthest neighbor). Distance between clusters
#'     is the maximum distance between any two points. Tends to produce compact clusters
#'     with small diameters.}
#'   \item{`"average"`}{Average linkage (UPGMA). Distance between clusters is the average
#'     of all pairwise distances. Compromise between single and complete linkage.}
#'   \item{`"mcquitty"`}{McQuitty's method (WPGMA). Similar to average linkage but uses
#'     weighted averages.}
#'   \item{`"median"`}{Median linkage (WPGMC). Uses median of pairwise distances.}
#'   \item{`"centroid"`}{Centroid linkage (UPGMC). Distance between cluster centroids.
#'     Can produce inversions in the dendrogram.}
#' }
#'
#' @section Performance and Optimization:
#' **Computational Complexity:**
#' \itemize{
#'   \item Distance matrix computation: O(p²n) where p = variables, n = observations
#'   \item Hierarchical clustering: O(p²log p)
#'   \item Total time complexity: O(p²n + p²log p)
#' }
#'
#'
#' **Caching:**
#' \itemize{
#'   \item `summary()` results are cached after first call
#'   \item Cache is invalidated when `fit()` or `cut_tree()` is called
#'   \item Repeated calls to `summary()` are instantaneous
#' }
#'
#' @section Data Requirements:
#' **Input Data Format:**
#' \itemize{
#'   \item Must be a `data.frame` or `matrix`
#'   \item Rows represent observations (samples)
#'   \item Columns represent variables to cluster
#'   \item Qualitative variables must be encoded as `factor`
#'   \item Quantitative variables must be `numeric`
#'   \item Missing values (NA) are handled with pairwise complete observations
#' }
#'
#' **Data Preprocessing:**
#' \itemize{
#'   \item Standardization is handled automatically for PCA
#'   \item No need to manually scale quantitative variables
#'   \item Ensure factors have meaningful levels
#'   \item Consider removing variables with >50% missing values
#'   \item One-hot encoded variables (0/1) are automatically detected
#' }
#'
#' @section Mathematical Details:
#' **Distance Computation:**
#'
#' For quantitative variables with metric "rsquare":
#' \deqn{d(X_i, X_j) = \sqrt{1 - cor(X_i, X_j)^2}}
#'
#' For quantitative variables with metric "r":
#' \deqn{d(X_i, X_j) = \sqrt{1 - cor(X_i, X_j)}}
#'
#' For qualitative variables:
#' \deqn{d(X_i, X_j) = 1 - V(X_i, X_j)}
#' where V is Cramer's V statistic.
#'
#' **Cramer's V:**
#' \deqn{V = \sqrt{\frac{\chi^2}{n \cdot min(r-1, c-1)}}}
#' where n = sample size, r = number of rows, c = number of columns in contingency table.
#'
#' **Correlation Ratio (η²):**
#' \deqn{\eta^2 = \frac{SSB}{SST} = \frac{\sum_{k} n_k(\bar{y}_k - \bar{y})^2}{\sum_{i} (y_i - \bar{y})^2}}
#' Used to measure association between categorical and continuous variables.
#'
#' @section Output Components:
#' **Summary Output:**
#'
#' `clust_summary` data frame contains:
#' \describe{
#'   \item{`cluster`}{Cluster identifier (integer)}
#'   \item{`n_members`}{Number of variables in the cluster}
#'   \item{`var_explained`}{Eigenvalue of first principal component (variance explained)}
#'   \item{`prop_explained`}{Proportion of variance explained per variable}
#' }
#'
#' `clust_members` data frame contains:
#' \describe{
#'   \item{`cluster`}{Assigned cluster number}
#'   \item{`own_cluster_R2`}{Squared correlation/association with own cluster centroid}
#'   \item{`next_closest_R2`}{Squared correlation/association with nearest alternative cluster}
#'   \item{`1 - R2_ratio`}{Quality metric: (1 - own_R2) / (1 - next_R2). Lower is better.}
#' }
#'
#' `centroids_correlations` data frame contains:
#' \describe{
#'   \item{`cluster A`}{First cluster in comparison}
#'   \item{`cluster B`}{Second cluster in comparison}
#'   \item{`correlation`}{Pearson correlation between cluster centroids}
#'   \item{`squared_correlations`}{R² between cluster centroids}
#' }
#'
#'
#' @examples
#' \dontrun{
#' library(datasets)
#' data(mtcars)
#'
#' hc <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
#' hc$fit(mtcars)
#'
#' # Visualize dendrogram
#' hc$plot_dendrogram()
#'
#' # Find optimal number of clusters
#' hc$plot_agg_levels()
#'
#' # Cut tree into 3 clusters
#' hc$cut_tree(k = 3)
#'
#' # Analyze results
#' results <- hc$summary()
#' print(results$clust_summary)
#' print(results$clust_members)
#'
#' # Quality assessment
#' hc$plot_silhouette()
#' hc$mds_projection()
#' }
#'
#'
#'
#' @references
#' Vigneau, E., & Qannari, E. M. (2003).
#' Clustering of variables around latent components.
#' Communications in Statistics - Simulation and Computation, 32(4), 1131-1150.
#' \doi{10.1081/SAC-120023882}
#'
#' @family clustering classes
#'
#' @export
HClustVar <- R6::R6Class(
  "HClustVar",
  inherit = ClusteringBase,

  # ==========================================================================
  # PRIVATE FIELDS
  # ==========================================================================
  private = list(

    # ==== Parameters

    # Metric used to compute distances between variables. (correlation or square correlation).
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

    # Store the result of the summary method.
    .summary_results = NULL,


    # ==========================================================================
    # PRIVATE METHODS
    # ==========================================================================

    # -----------------------------------------------------------------------
    # Check input method
    # -----------------------------------------------------------------------


    #' @description
    #' Checks that vartype, dist.metric, and cah.method have valid values
    #' and are compatible with each other. Called by initialize().
    #'
    #' @param vartype Character. Variable type ("quant", "qual", "mixed", "auto")
    #' @param dist.metric Character or NULL. Distance metric ("r", "rsquare", or NULL)
    #' @param cah.method Character. Clustering method (e.g., "ward.D")
    #'
    #' @return None. Stops with error if validation fails, warns if parameters incompatible.
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

    #' @description
    #' Automatically detect variable type from data
    #'
    #' Determines whether the dataset contains quantitative, qualitative, or
    #' mixed variables, and adjusts the distance metric accordingly.
    #'
    #' Called by \code{fit()} when \code{vartype = "auto"} to automatically
    #' determine the appropriate clustering method.
    #'
    #' **Detection logic:**
    #' \itemize{
    #'   \item No qualitative variables: "quant"
    #'   \item No quantitative variables: "qual"
    #'   \item Both types present: "mixed"
    #' }
    #'
    #' **Side effects:**
    #' \itemize{
    #'   \item Sets \code{private$.dist.metric = "rsquare"} if NULL for quantitative data
    #'   \item Warns and corrects invalid metric for quantitative data
    #'   \item Warns if metric specified for qualitative data (will be ignored)
    #' }
    #'
    #' @return Character string: "quant", "qual", or "mixed"
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

    #' @description
    #' Compute dissimilarity matrix for hierarchical clustering
    #' Calculates a distance matrix appropriate for the variable type
    #' (quantitative, qualitative, or mixed) to be used in hierarchical clustering.
    #'
    #' **For quantitative variables (vartype = "quant"):**
    #' \itemize{
    #'   \item Computes correlation matrix between variables
    #'   \item Applies metric transformation (r or r²)
    #'   \item Converts to dissimilarity: \code{dist = sqrt(1 - cor)}
    #' }
    #'
    #' **For qualitative variables (vartype = "qual"):**
    #' \itemize{
    #'   \item Computes Cramer's V matrix (association between categorical variables)
    #'   \item Converts to dissimilarity: \code{dist = 1 - cramer_v}
    #' }
    #'
    #' **For mixed variables (vartype = "mixed"):**
    #' \itemize{
    #'   \item Discretizes quantitative variables into 4 quantiles
    #'   \item Treats all variables as categorical
    #'   \item Computes Cramer's V matrix
    #'   \item Converts to dissimilarity: \code{dist = 1 - cramer_v}
    #' }
    #'
    #' @return A \code{dist} object containing pairwise dissimilarities between
    #'   variables, suitable for use with \code{hclust()}.
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

    #' @description
    #' Discretize quantitative variables into quantile-based categories
    #'
    #' @description
    #' Converts numeric variables into ordered factors using quantile-based binning.
    #' One-hot encoded variables (0/1) are converted to factors without binning.
    #'
    #' @param df Data frame containing variables to discretize
    #' @param quanti_index Integer vector of column indices for quantitative variables
    #' @param n_groups Number of quantile groups (bins) to create
    #'
    #' @return Data frame with specified columns converted to factors (labeled Q1, Q2, ...)
    #'
    #' Used internally for mixed variable clustering to treat all variables as categorical.
    #' Handles edge cases: empty indices, one-hot encoding, duplicate quantiles (small datasets).
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

    #' @description
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

    #' @description
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


    #' @description
    #' Compute cluster centroids using factorial analysis
    #'
    #' Computes the centroid (synthetic variable) for each cluster by performing
    #' factorial analysis (PCA, MCA, or FAMD) and extracting the first principal
    #' component. Also stores the eigenvalues representing variance explained.
    #'
    #' **Process for each cluster:**
    #' \enumerate{
    #'   \item Extracts variables belonging to the cluster
    #'   \item Performs appropriate factorial analysis:
    #'     \itemize{
    #'       \item PCA: All numeric variables (standardized)
    #'       \item MCA: All categorical variables
    #'       \item FAMD: Mixed numeric and categorical variables
    #'     }
    #'   \item Extracts first principal component as centroid
    #'   \item Stores eigenvalue (variance explained by centroid)
    #' }
    #'
    #' **Orientation correction (PCA only):**
    #'
    #' When using correlation distance (dist.metric = "r"), ensures the centroid
    #' is positively correlated with cluster variables on average by inverting
    #' the sign if mean correlation is negative.
    #'
    #' @return None (invisible NULL). Results stored in private fields.
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

        # Extract first component
        if (is.matrix(res$ind$coord)) {
          first_component <- res$ind$coord[, 1]
        } else if (is.data.frame(res$ind$coord)) {
          first_component <- res$ind$coord[[1]]
        } else {
          first_component <- as.numeric(res$ind$coord)
        }

        # ═══════════════════════════════════════════════════════════════
        # Correct orientation if needed
        # ═══════════════════════════════════════════════════════════════

        if (all(sapply(sub_data, is.numeric)) &&
            !is.null(self$dist.metric) &&
            !is.na(self$dist.metric) &&
            self$dist.metric == "r") {

            # Check correlation with variables.
            correlations <- apply(sub_data, MARGIN = 2, function(x) {
              cor(x, first_component, use = "pairwise.complete.obs")
            })


            mean_correlation <- mean(correlations, na.rm = TRUE)

            # If mean correlation if < 0, correct the sign of the latent component.
            if (!is.na(mean_correlation) && mean_correlation < 0) {
              first_component <- -first_component
            }
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

    # -----------------------------------------------------------------------
    # Correlation ratio
    # -----------------------------------------------------------------------

    #' @description
    #' Compute the correlation ratio (eta-squared) between categorical and numeric variables
    #'
    #' Calculates the correlation ratio (η²), a measure of association between
    #' a qualitative variable and a quantitative variable. This is the categorical
    #' equivalent of R² for continuous variables.
    #'
    #' The correlation ratio measures the proportion of variance in the quantitative
    #' variable that is explained by the categorical variable:
    #'
    #' \deqn{\eta^2 = \frac{SSB}{SST}}
    #'
    #' where:
    #' \itemize{
    #'   \item SSB: Between-group sum of squares
    #'   \item SST: Total sum of squares
    #' }
    #'
    #' **Usage Context:**
    #' Used internally for:
    #' \itemize{
    #'   \item Computing similarity between qualitative and quantitative variables
    #'   \item Variable assignment in mixed-type clustering
    #'   \item Summary statistics calculation for mixed clusters
    #' }
    #'
    #' **Interpretation:**
    #' \itemize{
    #'   \item η² = 0: No association (category doesn't explain the numeric variable)
    #'   \item η² = 1: Perfect association (category fully explains the numeric variable)
    #'   \item 0 < η² < 1: Partial association (higher values = stronger relationship)
    #' }
    #'
    #' @param quali A vector (will be coerced to factor) representing the categorical variable
    #' @param quanti A vector (will be coerced to numeric) representing the quantitative variable
    #'
    #' @return Numeric value between 0 and 1 representing the correlation ratio (η²)
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
    },

    # -----------------------------------------------------------------------
    # Silhouette computation method
    # -----------------------------------------------------------------------

    #' @description
    #' Compute silhouette coefficients for variable clustering
    #'
    #' Calculates silhouette coefficients for each variable based on their
    #' R² with own cluster and nearest cluster (from summary data).
    #'
    #' The silhouette coefficient is calculated as:
    #' s(i) = (own_cluster_R2 - next_closest_R2) / max(own_cluster_R2, next_closest_R2)
    #'
    #' The silhouette coefficient ranges from -1 to 1:
    #' - Values close to 1: variable is well matched to its cluster
    #' - Values close to 0: variable is on the border between two clusters
    #' - Negative values: variable might be assigned to the wrong cluster
    #'
    #' @return A numeric vector of silhouette coefficients, one per variable
    compute_silhouette = function() {

      # Check prerequisites
      if (!self$fitted) {
        stop("Model must be fitted before computing silhouette. Use fit() first.")
      }

      if (is.null(self$labels)) {
        stop("Tree must be cut before computing silhouette. Use cut_tree() first.")
      }

      # If only one cluster, silhouette is undefined
      if (self$n_clusters == 1) {
        stop("Silhouette is undefined for a single cluster")
      }

      # Get summary data (contains own_cluster_R2 and next_closest_R2)
      summary_data <- self$summary()
      clust_members <- summary_data$clust_members

      # Calculate silhouette coefficient
      # s(i) = (a(i) - b(i)) / max(a(i), b(i))
      # where a(i) = own_cluster_R2 and b(i) = next_closest_R2
      silhouette_coef <- with(clust_members, {
        (own_cluster_R2 - next_closest_R2) / pmax(own_cluster_R2, next_closest_R2)
      })

      # Handle edge cases (when both R2 are 0)
      silhouette_coef[is.nan(silhouette_coef)] <- 0
      silhouette_coef[is.infinite(silhouette_coef)] <- 0

      # Set names
      names(silhouette_coef) <- rownames(clust_members)

      return(silhouette_coef)
    }


  ),

  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================

  public = list(

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------

    #' @description
    #' Initialize a new class instance
    #'
    #' This method initializes the main parameters of the object, including the variable type
    #' (`vartype`), distance metric (`dist.metric`), and hierarchical clustering method
    #' (`cah.method`).
    #'
    #' If `vartype = "auto"`, automatic variable selection is activated based
    #' on the data frame properties that will be used for fit.
    #' If `vartype = "quant"` and `dist.metric` is `NULL`, the default metric `"rsquare"` is used.
    #'
    #' @param vartype Character string specifying the type of variables. Options: "auto" (enables automatic variable selection), "quant" (quantitative variables), "qual" (qualitative variables), or "mixed" (both types). Default is "auto".
    #' @param dist.metric Character string or NULL specifying the distance metric. If NULL and vartype = "quant", default is "rsquare". Options: "r" (correlation) or "rsquare" (squared correlation). Default is NULL.
    #' @param cah.method Character string specifying the hierarchical clustering method passed to \code{hclust}. Default is "ward.D".
    #'
    #'
    #' @details
    #' This method updates the following private fields:
    #' \itemize{
    #'   \item `private$.auto_var_selection`: flag indicating automatic variable selection.
    #'   \item `private$.dist.metric`: selected distance metric.
    #'   \item `private$.vartype`: type of variables.
    #'   \item `private$.cah.method`: hierarchical clustering method.
    #' }
    #'
    #' It also performs parameter validation through the private method
    #' `private$check_input()`.
    #'
    #' @return None
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

    #' @description
    #' Fit hierarchical clustering on variables
    #'
    #' Computes the hierarchical clustering based on the variable type (quantitative,
    #' qualitative, or mixed) and the specified distance metric.
    #'
    #' @param data data.frame or matrix containing the variables to cluster. Be careful to set your qualitatives or dummified data as factors.
    #'
    #' - If `vartype = "auto"`, the method automatically selects 'quant', 'qual', or 'mixed'
    #'   based on the input data.
    #' - For quantitative variables, computes correlation or squared correlation as dissimilarity.
    #' - For qualitative or mixed variables, transforms quantitative columns into quantiles (discretization)
    #'   and computes Cramer's V dissimilarity matrix.
    #' - Finally, applies hierarchical clustering (`hclust`) using the method defined in `.cah.method`.
    #' - Reset the summary cache.
    #'
    #' @return None. The clustering tree is stored internally in `private$.tree`.
    fit = function(data) {

      # reset the summary cache.
      private$.summary_results <- NULL

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

    #' @description
    #' Cut the hierarchical clustering tree into clusters
    #'
    #' Cuts the hierarchical clustering tree of a fitted model into clusters and compute their centroids.
    #' This method updates the `labels` and `n_clusters`attributes of the object with cluster assignments.
    #'
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
    #' This method reset the summary() cache.
    #'
    #' @param k Integer, optional. Number of clusters to generate.
    #' @param h Numeric, optional. Height at which to cut the tree.
    #'
    #' @return A vector of cluster labels corresponding to each observation. The labels are also
    #'         stored in the object's `labels` attribute.
    #'
    #'
    cut_tree = function(k = NULL, h = NULL) {

      # Check if the model is fitted.
      if (!self$fitted) {
        stop("Your model should be fitted on data first.")
      }

      # resset the summary cache.
      private$.summary_results <- NULL

      # Update the labels attribute.
      self$labels <- cutree(self$tree, k, h)

      # update the number of clusters
      self$n_clusters <- length(unique(self$labels))

      # Compute the centroids of the clusters.
      private$compute_centroids()

      # reinitialise the summary cache.
      private$.summary_results <- NULL

      return(self$labels)

    },

    # -----------------------------------------------------------------------
    # plot dendrogram method
    # -----------------------------------------------------------------------

    #' @description
    #' Plot dendrogram of the hierarchical clustering
    #'
    #' Displays the dendrogram of the computed hierarchical clustering tree,
    #' with optional colored rectangles highlighting the cluster groups.
    #'
    #' The dendrogram visualizes the hierarchical relationships between variables
    #' obtained from the hierarchical clustering algorithm.
    #' If a number of clusters `k` is provided, the function colors and outlines
    #' the corresponding groups in the dendrogram using distinct colors.
    #'
    #' @param k (integer) Optional. Number of clusters to highlight in the dendrogram.
    #'        If `k = 0` (default), the dendrogram is drawn without colored groups.
    #'
    #' @return None. Generates a plot.
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

    #' @description
    #' Predict cluster membership for new variables
    #'
    #' Assigns new variables to existing clusters based on their similarity
    #' with the variables used to build the clustering model.
    #'
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
    predict = function(new_data) {

      # ============
      # 1. Prerequisites validation
      # ============

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

      # ========
      # 2. INITIALIZATION
      # ========

      result <- integer(ncol(new_data))
      names(result) <- colnames(new_data)


      # If there's only one cluster, all variables belong to it, then return.
      if (self$n_clusters == 1) {
        result[] <- 1
        return(result)
      }

      # =========
      # 3. PREDICTION DEPENDING OF CAH METHOD
      # =========

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
              "mcquitty" = mean(cluster_similarities),     # UPGMC (Unweighted mean) -> Same as average.
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
    # Aggregation level plot
    # -----------------------------------------------------------------------

    #' @description
    #' Plot aggregation levels to determine optimal number of clusters
    #'
    #' Creates a diagnostic plot showing the hierarchical clustering aggregation
    #' heights as a function of the number of clusters. This plot helps determine
    #' the optimal number of clusters using the "elbow method" based on
    #' perpendicular distance.
    #'
    #' **Elbow Detection:**
    #'
    #' The optimal number of clusters is detected using perpendicular distance:
    #' \deqn{d_i = \frac{|(y_2 - y_1)x_i - (x_2 - x_1)y_i + x_2y_1 - y_2x_1|}{\sqrt{(y_2 - y_1)^2 + (x_2 - x_1)^2}}}
    #'
    #' where:
    #' \itemize{
    #'   \item \code{(x_i, y_i)}: Normalized coordinates of point i
    #'   \item \code{(x_1, y_1)}: First point (minimum clusters)
    #'   \item \code{(x_2, y_2)}: Last point (maximum clusters)
    #'   \item \code{k_optimal}: Point with maximum perpendicular distance
    #' }
    #'
    #' **Normalization:**
    #'
    #' Both axes are normalized to ensure scale-invariance:
    #' \itemize{
    #'   \item X-axis (cluster numbers) normalized to 0, 1
    #'   \item Y-axis (aggregation heights) normalized to 0, 1
    #'   \item This ensures equal weight for both dimensions in distance calculation
    #'   \item Results are independent of original measurement units
    #' }
    #'
    #' **Interpreting the Plot:**
    #'
    #' The ideal number of clusters is found at the "elbow" of the curve:
    #' \itemize{
    #'   \item Sharp initial decrease in aggregation height (distinct clusters)
    #'   \item Curve flattens after the elbow (forced merging of similar clusters)
    #'   \item Red vertical line marks the detected optimal k
    #'   \item Gray diagonal line shows the reference baseline
    #'   \item Look for the point furthest from the baseline
    #' }
    #'
    #' @param max_cluster Maximum number of clusters to display in the plot.
    #'   Default is 0 (all clusters). Must be between 0 and the maximum possible
    #'   number of clusters (number of observations). If 0, displays all available
    #'   cluster levels. If specified, limits the plot to the first \code{max_cluster}
    #'   levels for better visualization of early clustering structure.
    #'
    #' @return None. Produces a plot as a side effect showing:
    #'   \itemize{
    #'     \item Aggregation heights curve with points and lines
    #'     \item Gray diagonal reference line connecting endpoints
    #'     \item Red vertical line at optimal k
    #'     \item Red point highlighting the elbow position
    #'     \item Legend displaying the optimal k value
    #'   }
    plot_agg_levels = function(max_cluster = 0) {

      hc <- private$.tree
      # Get the height values
      heights <- rev(hc$height)
      n_clusters <- 1:(length(heights) + 1)
      n_clusters_plot <- NULL

      # Correction of the cluster numbers
      if (max_cluster < 0 || max_cluster > length(n_clusters)) {
        stop(paste0("max_cluster can't be < 0 or > of the maximum clusters possible: ", length(n_clusters)))

        # Default value: use all clusters
      } else if (max_cluster == 0) {

        # Reverse the n_clusters (last values beeing the highets aggregation level)
        n_clusters_plot <- n_clusters[-length(n_clusters)]
        heights_plot <- heights

      } else {

        # Limit to max_cluster clusters
        n_clusters_plot <- n_clusters[1:max_cluster]
        heights_plot <- heights[1:max_cluster]
      }

      # Normalize data to have comparable scales
      x <- n_clusters_plot
      y <- heights_plot

      # Normalization between 0 and 1
      x_norm <- (x - min(x)) / (max(x) - min(x))
      y_norm <- (y - min(y)) / (max(y) - min(y))

      # Calculate perpendicular distance from each point to the line
      # connecting the first and last point
      x1 <- x_norm[1]
      y1 <- y_norm[1]
      x2 <- x_norm[length(x_norm)]
      y2 <- y_norm[length(y_norm)]

      distances <- sapply(1:length(x_norm), function(i) {
        # Distance from a point to a line
        abs((y2 - y1) * x_norm[i] - (x2 - x1) * y_norm[i] + x2 * y1 - y2 * x1) /
          sqrt((y2 - y1)^2 + (x2 - x1)^2)
      })

      # The elbow is the point with maximum distance
      k_optimal <- n_clusters_plot[which.max(distances)]

      # Plot
      plot(n_clusters_plot, heights_plot,
           type = "b", pch = 19,
           xlab = "Number of clusters",
           ylab = "Aggregation level",
           main = "Elbow method (perpendicular distance)",
           xaxt = "n")
      axis(1, at = n_clusters_plot)
      grid()

      # Line connecting the endpoints
      segments(n_clusters_plot[1], heights_plot[1],
               n_clusters_plot[length(n_clusters_plot)], heights_plot[length(heights_plot)],
               col = "gray", lty = 2)

      # Mark the optimal k
      abline(v = k_optimal, col = "red", lty = 2, lwd = 2)
      points(k_optimal, heights_plot[which(n_clusters_plot == k_optimal)],
             col = "red", pch = 19, cex = 2)

      legend("topright",
             legend = paste("Optimal k =", k_optimal),
             col = "red", lty = 2, lwd = 2, bty = "n")
    },

    # -----------------------------------------------------------------------
    # Silhouette plot method
    # -----------------------------------------------------------------------

    #' @description
    #' Plot silhouette diagram for variable clustering
    #'
    #' Creates a professional silhouette plot showing how well each variable fits
    #' within its assigned cluster.
    #'
    #' @param main Character. Title for the plot. Default: "Silhouette Analysis"
    #' @param colors Vector of colors for each cluster. Default: professional palette
    #' @param show_values Logical. Display silhouette values on bars? Default: TRUE
    #' @param cex.names Numeric. Size of variable names. Default: 0.7
    #' @param sort_desc Logical. Sort by silhouette descending within clusters? Default: FALSE
    #'
    #' @return Invisibly returns the silhouette coefficients
    plot_silhouette = function(main = "Silhouette Analysis",
                               colors = NULL,
                               show_values = TRUE,
                               cex.names = 0.7,
                               sort_desc = FALSE) {

      # Compute silhouette coefficients
      sil_coef <- private$compute_silhouette()

      # Get summary for cluster information
      summary_data <- self$summary()
      clust_members <- summary_data$clust_members

      # Colors
      if (is.null(colors)) {
        if (self$n_clusters <= 8) {

          # Colorbrewer Set2 palette
          colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                      "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")[1:self$n_clusters]
        } else {
          colors <- rainbow(self$n_clusters, s = 0.6, v = 0.8)
        }
      }

      # Prepare data for plotting
      sil_data <- data.frame(
        variable = rownames(clust_members),
        cluster = clust_members$cluster,
        silhouette = sil_coef,
        stringsAsFactors = FALSE
      )

      # Sort by cluster, then by silhouette
      if (sort_desc) {
        sil_data <- sil_data[order(-sil_data$cluster, -sil_data$silhouette), ]
      } else {
        sil_data <- sil_data[order(-sil_data$cluster, sil_data$silhouette), ]
      }
      sil_data$index <- 1:nrow(sil_data)

      # Calculate statistics
      avg_sil_per_cluster <- tapply(sil_data$silhouette, sil_data$cluster, mean)
      overall_avg_sil <- mean(sil_data$silhouette)

      # Setup plot layout with proper margins
      old_par <- par(
        mar = c(5, 11, 4, 2),
        bg = "white",
        family = "sans"
      )
      on.exit(par(old_par))

      # Calculate plot limits
      x_min <- min(-0.15, min(sil_data$silhouette) - 0.05)
      x_max <- max(1.05, max(sil_data$silhouette) + 0.05)

      # Create empty plot
      plot(
        NA,
        xlim = c(x_min, x_max),
        ylim = c(0, nrow(sil_data) + 2),
        xlab = "",
        ylab = "",
        main = "",
        yaxt = "n",
        xaxt = "n",
        bty = "n"
      )

      # Add custom title
      mtext(main, side = 3, line = 2.5, cex = 1.4, font = 2, col = "gray20")
      mtext(
        sprintf("Average silhouette width: %.3f", overall_avg_sil),
        side = 3,
        line = 1,
        cex = 1,
        col = if(overall_avg_sil > 0.5) "darkgreen" else if(overall_avg_sil > 0.25) "orange" else "red",
        font = 2
      )

      # Add custom x-axis
      axis_breaks <- seq(ceiling(x_min * 10) / 10, floor(x_max * 10) / 10, by = 0.2)
      axis(1,
           at = axis_breaks,
           las = 1,
           col = "gray50",
           col.axis = "gray30",
           cex.axis = 1)
      mtext("Silhouette Coefficient", side = 1, line = 3.5, cex = 1.1, col = "gray20", font = 2)

      # Add subtle grid
      abline(v = seq(ceiling(x_min * 10) / 10, floor(x_max * 10) / 10, by = 0.1),
             col = "gray90", lty = 1)

      # Add reference lines
      abline(v = 0, col = "gray40", lty = 2, lwd = 1.5)
      abline(v = overall_avg_sil, col = "#E41A1C", lty = 2, lwd = 2.5)

      # Add label for average line
      text(
        x = overall_avg_sil,
        y = nrow(sil_data) + 1,
        labels = "Average",
        col = "#E41A1C",
        cex = 0.8,
        font = 2,
        pos = 3
      )

      # Plot bars for each cluster
      for (k in 1:self$n_clusters) {
        cluster_data <- sil_data[sil_data$cluster == k, ]

        # Draw bars and labels
        for (i in 1:nrow(cluster_data)) {
          y_pos <- cluster_data$index[i]
          sil_val <- cluster_data$silhouette[i]

          # Determine bar color intensity based on value
          bar_color <- colors[k]
          if (sil_val < 0) {
            bar_color <- adjustcolor(bar_color, alpha.f = 0.4)
          } else if (sil_val < 0.25) {
            bar_color <- adjustcolor(bar_color, alpha.f = 0.6)
          }

          # Draw silhouette bar with border
          rect(
            xleft = 0,
            xright = sil_val,
            ybottom = y_pos - 0.38,
            ytop = y_pos + 0.38,
            col = bar_color,
            border = adjustcolor(colors[k], alpha.f = 0.9),
            lwd = 0.8
          )

          # Add variable name
          text(
            x = x_min + 0.01,
            y = y_pos,
            labels = cluster_data$variable[i],
            pos = 4,
            cex = cex.names,
            col = "gray20",
            font = if(sil_val > avg_sil_per_cluster[k]) 2 else 1
          )

          # Optionally add silhouette value
          if (show_values && abs(sil_val) > 0.08) {
            text(
              x = sil_val,
              y = y_pos,
              labels = sprintf("%.2f", sil_val),
              pos = if(sil_val > 0) 4 else 2,
              cex = 0.6,
              col = "gray40"
            )
          }
        }
      }

      # Add subtle footer
      mtext(
        sprintf("%d variables • %d clusters", nrow(sil_data), self$n_clusters),
        side = 1,
        line = -1.5,
        adj = 0.98,
        cex = 0.75,
        col = "gray50"
      )

      invisible(sil_coef)
    },

    # -----------------------------------------------------------------------
    # Multidimensional scaling
    # -----------------------------------------------------------------------


    #' @description
    #' Plot multidimensional scaling projection of variables
    #'
    #' Creates a 2D visualization of variable relationships using classical
    #' multidimensional scaling (MDS) on the distance matrix. Variables are
    #' projected onto the first two principal coordinates, with automatic
    #' cluster coloring.
    #'
    #' **Algorithm:**
    #'
    #' Uses classical multidimensional scaling (\code{cmdscale}) to project
    #' the dissimilarity matrix into a 2D space while preserving pairwise
    #' distances as much as possible.
    #'
    #' **Interpretation:**
    #' \itemize{
    #'   \item Variables close together have high similarity
    #'   \item Variables far apart have low similarity
    #'   \item Axes represent principal coordinates (not directly interpretable)
    #'   \item Percentage on each axis shows variance explained by that dimension
    #' }
    #'
    #' **Quality Assessment:**
    #' \itemize{
    #'   \item Cumulative variance explained (title): Total information retained in 2D
    #'   \item High percentage (>70%): Good 2D representation of distances
    #'   \item Low percentage (<50%): Distances may be distorted in 2D projection
    #' }
    #'
    #' **Visualization:**
    #' \itemize{
    #'   \item If tree has been cut: Variables colored by cluster assignment
    #'   \item If tree not cut: All variables in black
    #'   \item Variable names displayed at their coordinates
    #'   \item Aspect ratio = 1 (equal scaling on both axes)
    #' }
    #'
    #' @return None. Produces a plot as a side effect.
    mds_projection = function() {

      if (!self$fitted) {
        stop("Model has to be fitted to compute mds projection")
      }

      mds <- cmdscale(d = self$dist.matrix, k = 2, eig = TRUE)

      # Detect eigens values negativity.
      if (any(mds$eig < 0)) {
        warning("One of the eigen value may be negative, the % of information restitued may be altered.")
      }

      # If one of the 2 first eigens values are negative, stop the method.
      if (any(mds$eig[1:2] < 0)) {
        stop("The eigen value are negative, and the variables could not be projected on the component.")
      }

      # Check the % of information restitued.
      eigens <- mds$eig[mds$eig > 0]

      # Cumulative quality
      ceigens <- cumsum(eigens) / sum(eigens)
      quality <- eigens / sum(eigens)

      # Plot the representation
      plot(mds$points[,1],mds$points[,2],type="n",cex=0.3,asp=1, xlab = paste("1st component", round(quality[1] * 100, 2), "%"), ylab=paste("2nd component", round(quality[2] * 100, 2), "%"))

      if (!is.null(self$labels)) {

        # Create a color palette.
        cluster_levels <- unique(self$labels)
        colors <- setNames(rainbow(length(cluster_levels)), cluster_levels)

        # Get the color for each label.
        point_colors <- colors[self$labels]

        text(mds$points[,1], mds$points[,2], labels = labels(self$dist.matrix), col = point_colors)
      } else {
        text(mds$points[,1], mds$points[,2], labels = labels(self$dist.matrix))
      }

      title(main = paste("Multidimensional scaling - Projection", round(ceigens[2]*100,1), "%"))
    },



    # -----------------------------------------------------------------------
    # Summary method
    # -----------------------------------------------------------------------

    #' @description
    #' Produce a detailed summary of the variable clustering results
    #'
    #' Generates comprehensive statistics about the clustering results, including
    #' cluster-level summaries and individual variable membership details. Results
    #' are cached for performance optimization.
    #'
    #' This method produces two main data frames:
    #'
    #' **1. Cluster Summary (`clust_summary`)**
    #'
    #' Contains aggregate statistics for each cluster:
    #' \itemize{
    #'   \item \code{cluster}: Cluster identifier (1 to k)
    #'   \item \code{n_members}: Number of variables in the cluster
    #'   \item \code{var_explained}: Eigenvalue of the first principal component
    #'         (amount of variance explained by the cluster centroid)
    #'   \item \code{prop_explained}: Proportion of variance explained per variable
    #'         (var_explained / n_members)
    #' }
    #'
    #' **2. Cluster Members (`clust_members`)**
    #'
    #' Contains detailed information for each variable:
    #' \itemize{
    #'   \item \code{cluster}: Assigned cluster number
    #'   \item \code{own_cluster_R2}: Squared correlation (or correlation ratio)
    #'         between the variable and its own cluster centroid. Values close to 1
    #'         indicate strong membership.
    #'   \item \code{next_closest_R2}: Squared correlation (or correlation ratio)
    #'         with the nearest alternative cluster centroid
    #'   \item \code{1 - R2_ratio}: Quality metric calculated as
    #'         \code{(1 - own_cluster_R2) / (1 - next_closest_R2)}. Lower values
    #'         indicate better cluster assignment (variable is much closer to its
    #'         own cluster than to the next closest one).
    #' }
    #'
    #' **3. Correlation between latent component of cluster**
    #'
    #' Contains proximity between clusters:
    #' \itemize{
    #'   \item \code{cluster A}: First comparison cluster
    #'   \item \code{cluster B}: Second comparison cluster
    #'   \item \code{correlation}: Pearson correlation
    #'   \item \code{suqared-correlation}: R2
    #' }
    #'
    #' **Similarity Metrics Used:**
    #' \itemize{
    #'   \item Quantitative variables: Squared Pearson correlation (R²)
    #'   \item Qualitative variables: Correlation ratio (η²)
    #' }
    #'
    #' **Sorting and Display:**
    #' \itemize{
    #'   \item Variables are sorted by cluster, then by \code{own_cluster_R2} (descending)
    #'   \item Variables in bold (higher R²) are better representatives of their cluster
    #'   \item The \code{1 - R2_ratio} helps identify poorly classified variables
    #'         (high values suggest potential misclassification)
    #' }
    #'
    #' **Performance Optimization:**
    #' Results are cached in \code{private$.summary_results} after the first call.
    #' The cache is invalidated when:
    #' \itemize{
    #'   \item \code{fit()} is called (new data)
    #'   \item \code{cut_tree()} is called (different number of clusters)
    #' }
    #'
    #' If the model is not fitted or the tree has not been cut, prints an
    #' informative message and returns \code{invisible(self)}.
    #'
    #' **Cluster Quality:**
    #' \itemize{
    #'   \item High \code{prop_explained} (> 0.7): Homogeneous cluster
    #'   \item Low \code{prop_explained} (< 0.5): Heterogeneous cluster
    #' }
    #'
    #' **Variable Assignment Quality:**
    #' \itemize{
    #'   \item \code{own_cluster_R2 > 0.7}: Strong membership
    #'   \item \code{own_cluster_R2 < 0.5}: Weak membership
    #'   \item \code{1 - R2_ratio < 0.3}: Well classified
    #'   \item \code{1 - R2_ratio > 0.7}: Potentially misclassified
    #' }
    #'
    #' @return A list with three elements:
    #' \describe{
    #'   \item{\code{clust_summary}}{Data frame with cluster-level statistics}
    #'   \item{\code{clust_members}}{Data frame with variable-level statistics}
    #'   \item{\code{centroids_correlation}}{Data frame with correlation between
    #'   clusters}
    #' }
    summary = function() {

      # Return the result if already cached.
      if (!is.null(private$.summary_results)) {
        return(private$.summary_results)
      }

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

      # If only 1 cluster observable, do not compute summary
      if (self$n_clusters == 1) {
        warning("Cannot compute proximity with other clusters if only 1 available.")
        return(NULL)
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

      # Compute Cluster members table
      clust_members <- data.frame(
        cluster = self$labels[variables]
      )

      own_cluster_R2 <- numeric(length(variables))
      next_closest_R2 <- numeric(length(variables))

      # Iteration on all variables to calculate their properties.
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
          cor_clusters_temp <- apply(self$centroids, MARGIN = 2, FUN = function(centroid) {
            private$correlation_ratio(self$data[[variable]], centroid)
          })
        }

        own_cluster_R2[i] <- cor_clusters_temp[own_clust]
        next_closest_R2[i] <- max(cor_clusters_temp[-own_clust])
      } # end for

      # Update the dataframe.
      clust_members["own_cluster_R2"] <- round(own_cluster_R2, 2)
      clust_members["next_closest_R2"] <- round(next_closest_R2, 2)
      clust_members["1 - R2_ratio"] <- round((1 - clust_members["own_cluster_R2"]) /
                                               (1 - clust_members["next_closest_R2"]), 2)

      # ----------------
      # Silhouette scores based on R²
      # ----------------

      silhouette_values <- (own_cluster_R2 - next_closest_R2) /
        pmax(own_cluster_R2, next_closest_R2)

      clust_members["silhouette"] <- round(silhouette_values, 3)

      # Average silhouette score
      avg_silhouette <- mean(silhouette_values)

      # Silhouette by cluster
      cluster_silhouette <- sapply(1:self$n_clusters, function(k) {
        cluster_vars <- which(self$labels[variables] == k)
        if (length(cluster_vars) > 0) {
          mean(silhouette_values[cluster_vars])
        } else {
          NA
        }
      })

      clust_summary["avg_silhouette"] <- round(cluster_silhouette, 3)

      # Sort the clust_members df
      clust_members <- clust_members[order(clust_members$cluster, -clust_members$own_cluster_R2), ]

      # ----------------
      # Correlation between latent components
      # ----------------

      # Initialise the dataframe with all cluster pairs
      component_correlations <- as.data.frame(t(combn(colnames(self$centroids), 2)))
      colnames(component_correlations) <- c("cluster A", "cluster B")

      # iterate over all pairs to compute correlation on a "correlation column"
      component_correlations$correlation <- apply(component_correlations, 1, function(row) {
        round(cor(self$centroids[,row[1]], self$centroids[,row[2]], use = "complete.obs"), 3)
      })

      # Add a squared correlation column
      component_correlations$squared_correlations <- round(component_correlations$correlation**2, 3)

      # Prepare result object
      result <- list(
        clust_summary = clust_summary,
        clust_members = clust_members,
        centroids_correlations = component_correlations,
        avg_silhouette = round(avg_silhouette, 3),
        total_var_explained = round(percentage_var_explained * 100, 2)
      )

      # update the private attribute to cache the data.
      private$.summary_results <- result

      return(result)
    },

    # -----------------------------------------------------------------------
    # Print method
    # -----------------------------------------------------------------------


    #' @description
    #' Print method for HClustVar objects
    #'
    #' Displays a concise overview of the HClustVar object with key information.
    #'
    #' @return Invisibly returns the object itself.
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
      cat("Use $summary() for detailed information and cluster analysis\n")
      cat("Use $plot_dendrogram() to visualize the tree\n")
      cat("Use $plot_agg_levels() to visualize the inertia of clusters\n")
      cat("Use $plot_silhouette() to visualize the silhouette plot of clusters\n")
      cat("Use $mds_projection() to compute multidimensional scaling\n\n")

      invisible(self)
    }
  ),

  # ==========================================================================
  # ACTIVE BINDINGS (getters/setters)
  # ==========================================================================

  active = list(

    #' @field dist.metric Distance metric for quantitative variables.
    #'   Returns "r" (correlation), "rsquare" (squared correlation), or NULL
    #'   (for qualitative/mixed variables). Read-only.
    dist.metric = function() {return(private$.dist.metric)},

    #' @field vartype Type of variables in the dataset.
    #'   Returns "quant" (quantitative), "qual" (qualitative), "mixed", or "auto"
    #'   (before fitting). Read-only after auto-detection.
    vartype = function() {return(private$.vartype)},

    #' @field dist.matrix Dissimilarity matrix computed from the data.
    #'   Returns a \code{dist} object suitable for hierarchical clustering.
    #'   Available after \code{fit()} is called. Read-only.
    dist.matrix = function() {return(private$.dist.matrix)},

    #' @field tree Hierarchical clustering tree.
    #'   Returns an \code{hclust} object containing the dendrogram structure.
    #'   Available after \code{fit()} is called. Read-only.
    tree = function() {return(private$.tree)},

    #' @field cah.method Hierarchical clustering method used.
    #'   Returns one of: "ward.D", "ward.D2", "single", "complete", "average",
    #'   "mcquitty", "median", "centroid". Read-only.
    cah.method = function() {return(private$.cah.method)},

    #' @field summary_results Cached results from the \code{summary()} method.
    #'   Returns a list with \code{clust_summary} and \code{clust_members} data frames,
    #'   or NULL if summary has not been computed yet. Cache is invalidated when
    #'   \code{fit()} or \code{cut_tree()} is called. Read-only.
    summary_results = function() {return(private$.summary_results)},


    #' @field centroids Cluster centroids as synthetic variables.
    #'   Returns a matrix (n_observations × n_clusters) where each column represents
    #'   the first principal component of a cluster. Available after \code{cut_tree()}
    #'   is called. Can be set manually (advanced use only).
    #'
    #'   **Getter**: Returns the centroids matrix or NULL if tree not cut.
    #'
    #'   **Setter**: Assigns a custom centroids matrix. Must be:
    #'   \itemize{
    #'     \item A data.frame or matrix
    #'     \item Same number of rows as training data
    #'     \item Same number of columns as number of clusters
    #'   }
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

    #' @field clusters.eigen Eigenvalues (variance explained) for each cluster.
    #'   Returns a numeric vector where each element represents the variance explained
    #'   by the first principal component of the corresponding cluster. Values range
    #'   from 0 (no variance explained) to the number of variables in the cluster
    #'   (all variance explained). Available after \code{cut_tree()} is called.
    #'   Can be set manually (advanced use only).
    #'
    #'   **Getter**: Returns the eigenvalues vector or NULL if tree not cut.
    #'
    #'   **Setter**: Assigns custom eigenvalues. Must be:
    #'   \itemize{
    #'     \item A numeric vector
    #'     \item Length equal to the number of clusters
    #'   }
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

