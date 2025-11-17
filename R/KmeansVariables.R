#' K-means Variable Clustering
#'
#' @description
#' Implements a variable clustering algorithm using K-means adapted for variables.
#' This reallocation-based algorithm groups variables by maximizing their correlation
#' with cluster centroids (first principal component of each cluster).
#'
#' @details
#' The `KmeansVariables` class inherits from `ClusteringBase` and provides:
#' - Iterative reallocation of variables to clusters
#' - Distance calculation based on correlation with centroids (R² or signed r)
#' - Centroid computation via Principal Component Analysis (PCA)
#' - Multiple random initializations to avoid local minima
#' - Convergence monitoring and iteration tracking
#'
#' @section Algorithm:
#' The K-means for variables works as follows:
#' \enumerate{
#'   \item **Initialization**: Randomly assign variables to K clusters
#'   \item **Centroid Calculation**: For each cluster, compute the first principal component (via PCA) as the cluster centroid
#'   \item **Distance Calculation**: Compute distance between each variable and each centroid based on their correlation
#'   \item **Reallocation**: Reassign each variable to its nearest centroid
#'   \item **Convergence Check**: Repeat steps 2-4 until no reassignment occurs or max iterations reached
#' }
#'
#' @section Distance Metrics:
#' Two distance metrics are available for measuring dissimilarity between variables and centroids:
#'
#' **R²-based distance** (\code{distance_metric = "r_squared"}):
#' \deqn{d_{R^2}(X_j, c_k) = \sqrt{1 - \rho^2(X_j, c_k)}}
#' where \eqn{\rho(X_j, c_k)} is the Pearson correlation coefficient.
#' This metric groups variables with high absolute correlation (positive or negative).
#'
#' **Signed r-based distance** (\code{distance_metric = "r_signed"}):
#' \deqn{d_r(X_j, c_k) = \sqrt{1 - \rho(X_j, c_k)}}
#' This metric only groups variables with positive correlation, treating
#' negatively correlated variables as dissimilar.
#'
#' @section Objective Function:
#' The algorithm maximizes the within-cluster cohesion (inertia). For cluster \eqn{k},
#' the inertia is defined as the first eigenvalue \eqn{\lambda_k} from the PCA
#' of variables in \eqn{C_k}, which equals:
#' \deqn{I_k = \lambda_k = \sum_{X_j \in C_k} \rho^2(X_j, c_k)}
#' where \eqn{C_k} denotes the set of variables in cluster \eqn{k} and \eqn{c_k} is the
#' first principal component. The total inertia is \eqn{I = \sum_{k=1}^K I_k}.
#'
#' @section Convergence:
#' The algorithm stops when one of the following conditions is met:
#' \itemize{
#'   \item No variable changes cluster between two consecutive iterations
#'   \item The change in total inertia is less than \code{tol}
#'   \item Maximum number of iterations (\code{max_iter}) is reached
#' }
#'
#' @section Multiple Initializations:
#' To avoid local maxima, the algorithm runs \code{n_init} times with different
#' random initializations. The solution with the highest inertia is retained.
#' Set \code{random_state} to an integer for reproducible results.
#'
#' @examples
#' \dontrun{
#' # Example with synthetic data (3 clusters of correlated variables)
#' set.seed(123)
#' n <- 500
#'
#' # Create 3 groups of 5 correlated variables each
#' create_cluster <- function(n_vars, n_obs, cor_strength) {
#'   cov_matrix <- matrix(cor_strength, n_vars, n_vars)
#'   diag(cov_matrix) <- 1
#'   MASS::mvrnorm(n_obs, rep(0, n_vars), cov_matrix)
#' }
#'
#' X1 <- create_cluster(5, n, 0.7)
#' X2 <- create_cluster(5, n, 0.65)
#' X3 <- create_cluster(5, n, 0.6)
#'
#' data_test <- data.frame(X1, X2, X3)
#' colnames(data_test) <- paste0("Var", 1:15)
#'
#' # Option 1: Manual K selection (you choose K)
#' km <- KmeansVariables$new(n_clusters = 3, random_state = 42)
#' km$fit(data_test)
#'
#' # Display results
#' km$print()
#' km$summary()
#'
#' # Access properties
#' km$clusters      # Cluster assignments
#' km$inertia       # Total inertia
#' km$n_iter        # Number of iterations
#'
#' # Option 2: Automatic K selection (algorithm chooses optimal K)
#' km_auto <- KmeansVariables$new(n_clusters = "auto",
#'                                k_range = 2:8,
#'                                selection_method = "silhouette",
#'                                random_state = 42)
#' km_auto$fit(data_test)
#' # → Displays: "Optimal K selected: 3 (method: silhouette, avg_silhouette: 0.XXX)"
#'
#' # Option 3: Use signed correlation (separates positively/negatively correlated vars)
#' km_signed <- KmeansVariables$new(n_clusters = 3,
#'                                  distance_metric = "r_signed",
#'                                  random_state = 42)
#' km_signed$fit(data_test)
#' # Variables with negative correlation to centroid will be in different clusters
#'
#' # Calculate metrics
#' sil <- km$silhouette()                    # Calculate silhouette
#' ch <- km$calinski_harabasz()              # Calculate CH index
#' contrib <- km$contributions()             # Get variable contributions
#' cor_table <- km$correlation_table()       # Get correlation table
#'
#' # Visualization: 3 essential plots
#' km$plot_elbow(k_range = 2:8)              # Choose K (elbow method)
#' km$plot_projection(show_labels = TRUE)    # See variable groups (PCA circle)
#' km$plot_contributions(top_n = 5)          # Identify key variables per cluster
#'
#' # Predict on new variables
#' new_vars <- data.frame(
#'   NewVar1 = rnorm(n),
#'   NewVar2 = rnorm(n)
#' )
#' predictions <- km$predict(new_vars)
#' print(predictions$clusters)
#' 
#' # Example with real data (FIFA dataset)
#' # Assuming players_22_cleaned.csv is loaded
#' fifa_skills <- read.csv("players_22_cleaned.csv")
#' 
#' # Select skill variables only
#' skills <- fifa_skills[, grep("skill_|pace|shooting|passing", names(fifa_skills))]
#' 
#' # Cluster with 5 groups (attack, defense, technique, physical, goalkeeper)
#' km_fifa <- KmeansVariables$new(n_clusters = 5, random_state = 123)
#' km_fifa$fit(skills)
#' km_fifa$summary()
#' }
#'
#' @section Dependencies:
#' This class requires:
#' \itemize{
#'   \item The R6 package
#'   \item The parent class ClusteringBase
#'   \item Base R functions: prcomp(), cor(), sample()
#' }
#'
#' @param n_clusters Integer or "auto". Number of clusters. If "auto", optimal K selected during fit()
#' @param max_iter Integer. Maximum number of iterations. Default: 100
#' @param tol Numeric. Convergence tolerance. Default: 1e-4
#' @param n_init Integer. Number of random initializations. Default: 10
#' @param random_state Integer or NULL. Random seed for reproducibility. Default: NULL
#' @param k_range Integer vector. Range of K values to test when n_clusters = "auto". Default: 2:10
#' @param selection_method Character. Method for automatic K selection: "silhouette", "calinski", or "all". Default: "silhouette"
#' @param distance_metric Character. Metric for distance calculation:
#'        "r_squared" (based on R², groups variables with high absolute correlation) or
#'        "r_signed" (based on signed r, groups only positively correlated variables).
#'        Default: "r_squared"
#' @param data Data frame or matrix with variables to cluster
#' @param round_digits Integer. Number of decimal places for rounding. Default: 3
#' @param top_n Integer. Number of top variables to display. Default: 5
#' @param ... Additional arguments passed to plot functions
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new(n_clusters, max_iter, tol, n_init, random_state, k_range, selection_method, distance_metric)}}{
#'     Initialize a new KmeansVariables object. See constructor documentation for parameter details.
#'   }
#'   \item{\code{$fit(data)}}{
#'     Fit the K-means model to the data. Runs multiple random initializations and keeps best solution.
#'   }
#'   \item{\code{$predict(data)}}{
#'     Predict cluster assignments for new variables based on fitted centroids.
#'   }
#'   \item{\code{$print()}}{Display concise model information}
#'   \item{\code{$summary()}}{Display detailed model information with statistics}
#'   \item{\code{$silhouette()}}{Calculate silhouette coefficients for all variables}
#'   \item{\code{$calinski_harabasz()}}{Calculate Calinski-Harabasz index (cluster quality metric)}
#'   \item{\code{$intra_correlation()}}{Calculate mean intra-cluster correlations}
#'   \item{\code{$contributions()}}{Calculate variable contributions to each cluster}
#'   \item{\code{$correlation_table(round_digits = 3)}}{Get correlation matrix between variables and centroids}
#'   \item{\code{$plot_elbow(k_range = 2:10, ...)}}{Plot elbow curve for optimal K selection}
#'   \item{\code{$plot_contributions(top_n = 5, ...)}}{Plot top contributing variables per cluster}
#'   \item{\code{$plot_projection(...)}}{Plot PCA-style projection of variables and clusters}
#' }
#'
#' @field clusters Integer vector of cluster assignments for each variable
#' @field centroids Matrix of cluster centroids (PC1 scores for each cluster)
#' @field inertia Total inertia (sum of squared correlations with centroids)
#' @field n_iter Number of iterations performed in last fit
#' @field cluster_sizes Integer vector of number of variables per cluster
#' @field cluster_inertias Numeric vector of inertia per cluster
#'
#' @seealso
#' \code{\link{ClusteringBase}} for the parent class
#' \code{\link[stats]{prcomp}} for Principal Component Analysis
#' \code{\link[stats]{cor}} for correlation calculations
#'
#' @references
#' \itemize{
#'   \item Vigneau, E., & Qannari, E. M. (2003). Clustering of variables around
#'         latent components. Communications in Statistics-Simulation and Computation,
#'         32(4), 1131-1150.
#'   \item Chavent, M., Kuentz-Simonet, V., Labenne, A., & Saracco, J. (2012).
#'         ClustOfVar: An R Package for the Clustering of Variables.
#'         Journal of Statistical Software, 50(13), 1-16.
#' }
#'
#' @family clustering classes
#'
#' @export
KmeansVariables <- R6::R6Class(
  "KmeansVariables",
  inherit = ClusteringBase,
  
  # ==========================================================================
  # PRIVATE FIELDS
  # ==========================================================================
  private = list(
    
    # Algorithm parameters
    .max_iter = NULL,        # Maximum number of iterations
    .tol = NULL,             # Convergence tolerance
    .n_init = NULL,          # Number of random initializations
    .random_state = NULL,    # Seed for reproducibility
    .distance_metric = NULL,  # Metric for distance calculation: "r_squared" (R²) or "r_signed" (signed r)

    # Automatic K selection parameters
    .k_range = NULL,         # Range of K values to test for automatic selection
    .selection_method = NULL, # Method for automatic K selection ("silhouette", "calinski", "all")

    # Results storage
    .centroids = NULL,       # Matrix of cluster centroids (PCA components)
    .inertia = NULL,         # Total within-cluster inertia
    .n_iter = NULL,          # Number of iterations performed
    .cluster_inertias = NULL, # Inertia per cluster
    
    # ==========================================================================
    # PRIVATE METHODS - CORE ALGORITHM
    # ==========================================================================
    
    # Validate algorithm parameters
    #
    # Checks that max_iter, tol, n_init are valid positive values
    validate_params = function(max_iter, tol, n_init) {
      
      if (!is.numeric(max_iter) || max_iter < 1) {
        stop("'max_iter' must be a positive integer")
      }
      
      if (!is.numeric(tol) || tol <= 0) {
        stop("'tol' must be a positive number")
      }
      
      if (!is.numeric(n_init) || n_init < 1) {
        stop("'n_init' must be a positive integer")
      }
    },
    
    # Calculate centroid of a cluster using PCA
    #
    # Computes the first principal component of variables in a cluster.
    # This component serves as the cluster centroid.
    # Also returns PCA metrics (eigenvalue, correlation loadings) for inertia calculation.
    #
    # @param X Matrix or data.frame of observations × variables
    # @param var_indices Integer vector of column indices in the cluster
    # @return List with:
    #   - centroid: Numeric vector (PC1 scores) of length nrow(X)
    #   - eigenvalue: First eigenvalue (variance explained by PC1)
    #   - cor_loadings: Correlation loadings (for inertia calculation)
    #   - n_vars: Number of variables in cluster
    calculate_centroid = function(X, var_indices) {

      n_vars <- length(var_indices)

      # Handle single variable case
      if (n_vars == 1) {
        return(list(
          centroid = as.numeric(scale(X[, var_indices])),
          eigenvalue = 1,  # Single variable explains 100% of itself
          cor_loadings = 1,
          n_vars = 1
        ))
      }

      # Extract cluster variables
      X_cluster <- X[, var_indices, drop = FALSE]

      # Remove constant variables (zero variance)
      variances <- apply(as.matrix(X_cluster), 2, function(x) var(x))
      if (any(variances == 0)) {
        X_cluster <- X_cluster[, variances > 0, drop = FALSE]
      }

      # If only one variable left or no variables after filtering
      if (ncol(X_cluster) <= 1) {
        if (ncol(X_cluster) == 1) {
          return(list(
            centroid = as.numeric(scale(X_cluster[, 1])),
            eigenvalue = 1,
            cor_loadings = 1,
            n_vars = n_vars
          ))
        } else {
          # No variables left, return zeros
          return(list(
            centroid = rep(0, nrow(X_cluster)),
            eigenvalue = 0,
            cor_loadings = numeric(0),
            n_vars = n_vars
          ))
        }
      }

      # Perform PCA (scale variables first)
      pca_result <- tryCatch({
        prcomp(X_cluster, center = TRUE, scale. = TRUE)
      }, error = function(e) {
        stop(paste0("PCA failed for cluster ", k, ": ", e$message,
                    "\nThis may indicate numerical instability or insufficient data."))
      })

      # Extract first principal component (scores)
      centroid <- pca_result$x[, 1]

      # Extract eigenvalue (variance explained by PC1)
      eigenvalue <- pca_result$sdev[1]^2

      # Calculate correlation loadings: rotation * sqrt(eigenvalue)
      # These represent correlations between original variables and PC1
      cor_loadings <- pca_result$rotation[, 1] * pca_result$sdev[1]

      # Align centroid so that majority of variables are positively correlated
      # This is crucial for distance_metric = "r_signed" to work correctly
      if (sum(cor_loadings) < 0) {
        centroid <- -centroid
        cor_loadings <- -cor_loadings
      }

      return(list(
        centroid = as.numeric(centroid),
        eigenvalue = eigenvalue,
        cor_loadings = cor_loadings,  # Store with sign for proper alignment
        n_vars = n_vars
      ))
    },
    
    # Compute distances between variables and centroids
    #
    # Calculates the distance matrix based on correlation.
    # If distance_metric = "r_squared": Distance = sqrt(1 - cor(X_j, c_k)^2) (R²-based)
    # If distance_metric = "r_signed": Distance = sqrt(1 - cor(X_j, c_k)) (r-based, sign matters)
    #
    # @param X Matrix or data.frame (n observations × p variables)
    # @param centroids Matrix (n observations × K clusters)
    # @param return_correlations Logical. If TRUE, also return correlation matrix
    # @return Matrix (p variables × K clusters) of distances, or list with distances and correlations
    compute_distances = function(X, centroids, return_correlations = FALSE) {

      p <- ncol(X)
      K <- ncol(centroids)

      dist_matrix <- matrix(0, nrow = p, ncol = K)
      cor_matrix <- matrix(0, nrow = p, ncol = K)
      rownames(dist_matrix) <- rownames(cor_matrix) <- colnames(X)
      colnames(dist_matrix) <- colnames(cor_matrix) <- paste0("Cluster", 1:K)

      for (j in 1:p) {
        for (k in 1:K) {
          # Correlation between variable j and centroid k
          cor_val <- suppressWarnings(cor(X[, j], centroids[, k]))

          # Handle NA from zero variance: set distance to 1 (maximum)
          if (is.na(cor_val)) {
            cor_val <- 0
          }

          # Store correlation
          cor_matrix[j, k] <- cor_val

          # Calculate distance based on distance_metric
          if (private$.distance_metric == "r_squared") {
            # R²-based: groups variables with high |correlation| (positive OR negative)
            dist_matrix[j, k] <- sqrt(1 - cor_val^2)
          } else {
            # r-based: groups only variables with high positive correlation
            # Uses correlation sign, so negatively correlated variables are pushed away
            dist_matrix[j, k] <- sqrt(1 - cor_val)
          }
        }
      }

      if (return_correlations) {
        return(list(distances = dist_matrix, correlations = cor_matrix))
      } else {
        return(dist_matrix)
      }
    },
    
    # Assign variables to nearest centroids
    # 
    # @param dist_matrix Matrix (p variables × K clusters) of distances
    # @return Integer vector of cluster assignments (1 to K)
    assign_clusters = function(dist_matrix) {
      
      # For each variable, find cluster with minimum distance
      clusters <- apply(dist_matrix, 1, which.min)
      
      return(as.integer(clusters))
    },
    
    # Single run of K-means algorithm
    #
    # @param X Matrix or data.frame (n × p)
    # @return List with clusters, centroids, inertia, n_iter
    kmeans_single_run = function(X) {

      n <- nrow(X)
      p <- ncol(X)
      K <- private$.n_clusters

      # Random initialization
      clusters <- sample(1:K, p, replace = TRUE)

      # Ensure all clusters have at least one variable
      for (k in 1:K) {
        if (sum(clusters == k) == 0) {
          clusters[sample(1:p, 1)] <- k
        }
      }

      inertia_old <- -Inf

      for (iter in 1:private$.max_iter) {

        # Step 1: Calculate centroids and collect PCA metrics
        centroids <- matrix(0, nrow = n, ncol = K)
        colnames(centroids) <- paste0("Cluster", 1:K)
        pca_metrics <- vector("list", K)  # Store PCA results for each cluster

        for (k in 1:K) {
          var_indices <- which(clusters == k)
          pca_result <- private$calculate_centroid(X, var_indices)
          centroids[, k] <- pca_result$centroid
          pca_metrics[[k]] <- pca_result
        }

        # Step 2: Compute distances
        dist_matrix <- private$compute_distances(X, centroids)

        # Step 3: Reassign clusters
        clusters_new <- private$assign_clusters(dist_matrix)

        # Ensure no empty clusters
        for (k in 1:K) {
          if (sum(clusters_new == k, na.rm = TRUE) == 0) {
            # Assign farthest variable to empty cluster
            farthest_var <- which.max(apply(dist_matrix, 1, min))
            clusters_new[farthest_var] <- k
          }
        }

        # Step 4: Calculate total inertia
        inertia_new <- 0
        for (k in 1:K) {
          n_vars <- pca_metrics[[k]]$n_vars

          if (n_vars == 0) {
            next
          }

          # Cluster inertia = eigenvalue of first principal component
          inertia_new <- inertia_new + pca_metrics[[k]]$eigenvalue
        }

        # Check convergence
        if (isTRUE(all.equal(clusters_new, clusters))) {
          clusters <- clusters_new
          break
        }

        if (abs(inertia_old - inertia_new) < private$.tol) {
          clusters <- clusters_new
          break
        }

        clusters <- clusters_new
        inertia_old <- inertia_new
      }

      # Calculate final metrics using PCA results
      # Recalculate centroids and PCA metrics for final clusters
      centroids_final <- matrix(0, nrow = n, ncol = K)
      colnames(centroids_final) <- paste0("Cluster", 1:K)
      cluster_inertias <- numeric(K)
      names(cluster_inertias) <- paste0("Cluster", 1:K)

      final_inertia <- 0
      for (k in 1:K) {
        var_indices <- which(clusters == k)
        pca_result <- private$calculate_centroid(X, var_indices)
        centroids_final[, k] <- pca_result$centroid

        n_vars <- pca_result$n_vars
        cor_loadings <- pca_result$cor_loadings

        # Calculate cluster inertia
        cluster_inertias[k] <- pca_result$eigenvalue

        final_inertia <- final_inertia + cluster_inertias[k]
      }

      list(
        clusters = clusters,
        centroids = centroids_final,
        inertia = final_inertia,
        n_iter = iter,
        cluster_inertias = cluster_inertias
      )
    }
  ),
  
  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================
  
  public = list(
    
    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------
    
    # @title Constructor for KmeansVariables
    #
    # @description
    # Initializes a new instance of the `KmeansVariables` class.
    #
    # @param n_clusters Integer or "auto". Number of clusters.
    #        If "auto", the optimal K will be selected automatically during fit(). Default: 3.
    # @param max_iter Integer. Maximum number of iterations. Default: 100.
    # @param tol Numeric. Convergence tolerance. Default: 1e-4.
    # @param n_init Integer. Number of random initializations. Default: 10.
    # @param random_state Integer or NULL. Random seed for reproducibility. Default: NULL.
    # @param k_range Integer vector. Range of K values to test when n_clusters = "auto".
    #        Default: 2:10. Ignored if n_clusters is an integer.
    # @param selection_method Character. Method for automatic K selection:
    #        "silhouette" (maximize avg silhouette),
    #        "calinski" (maximize Calinski-Harabasz index),
    #        or "all" (combine both, prefer silhouette if disagreement).
    #        Default: "silhouette". Ignored if n_clusters is an integer.
    # @param distance_metric Character. Metric for distance calculation:
    #        "r_squared" - distance based on R² = sqrt(1 - cor²), groups variables
    #                      with high absolute correlation (positive or negative)
    #        "r_signed"  - distance based on signed r = sqrt(1 - cor), groups only
    #                      positively correlated variables
    #        Default: "r_squared"
    #
    # @return A new instance of `KmeansVariables`.
    initialize = function(n_clusters = 3, max_iter = 100, tol = 1e-4,
                          n_init = 10, random_state = NULL,
                          k_range = 2:10, selection_method = "silhouette",
                          distance_metric = "r_squared") {

      # Validate parameters
      private$validate_params(max_iter, tol, n_init)

      # Validate distance_metric
      if (!distance_metric %in% c("r_squared", "r_signed")) {
        stop("'distance_metric' must be either 'r_squared' or 'r_signed'")
      }
      private$.distance_metric <- distance_metric

      # Handle n_clusters = "auto"
      if (is.character(n_clusters) && n_clusters == "auto") {
        # Store "auto" as special value (will be resolved during fit)
        private$.n_clusters <- "auto"

        # Validate k_range
        if (!is.numeric(k_range) || length(k_range) < 2) {
          stop("'k_range' must be a numeric vector with at least 2 values")
        }
        if (any(k_range < 2)) {
          stop("All values in 'k_range' must be >= 2")
        }

        # Validate selection_method
        if (!selection_method %in% c("silhouette", "calinski", "all")) {
          stop("'selection_method' must be 'silhouette', 'calinski', or 'all'")
        }

        # Store auto-selection parameters
        private$.k_range <- as.integer(k_range)
        private$.selection_method <- selection_method

      } else {
        # Set n_clusters via parent class setter (includes validation)
        self$n_clusters <- n_clusters

        # No need for k_range/selection_method
        private$.k_range <- NULL
        private$.selection_method <- NULL
      }

      # Store algorithm parameters
      private$.max_iter <- as.integer(max_iter)
      private$.tol <- tol
      private$.n_init <- as.integer(n_init)
      private$.random_state <- random_state

      # Set random seed if provided
      if (!is.null(random_state)) {
        set.seed(random_state)
      }
    },
    
    # -----------------------------------------------------------------------
    # Fit method
    # -----------------------------------------------------------------------
    
    # Fit K-means clustering on variables
    #
    # @description
    # Performs K-means clustering on quantitative variables.
    # Runs multiple initializations and keeps the best solution.
    #
    # @param data data.frame or matrix containing quantitative variables to cluster.
    #
    # @details
    # - Validates that data contains at least 2 quantitative variables
    # - If the dataset contains qualitative variables, they are automatically ignored
    #   (a message is displayed to inform the user)
    # - Runs the algorithm `n_init` times with different random starts
    # - Keeps the solution with the lowest inertia
    # - Stores results in private fields accessible via active bindings
    #
    # @return None. Results are stored internally and accessible via getters.
    fit = function(data) {
      
      # Validate and load data
      self$load_and_check_data(data)
      
      # Ensure we have at least 2 quantitative variables
      self$validate_algorithm_requirements("quant")
      
      # Extract quantitative data only
      X_quanti <- self$get_quanti_data()

      # Inform user if qualitative variables were detected and ignored
      n_quali <- length(private$.quali_indices)
      n_quanti <- length(private$.quanti_indices)

      if (n_quali > 0) {
        message(sprintf(
          "Note: %d qualitative variable%s detected and ignored. Using %d quantitative variable%s for clustering.",
          n_quali,
          ifelse(n_quali > 1, "s", ""),
          n_quanti,
          ifelse(n_quanti > 1, "s", "")
        ))
      }

      # Update .data to contain only quantitative variables
      # This ensures labels length matches the data used for clustering
      private$.data <- X_quanti
      private$.quali_indices <- integer(0)  # No more quali vars

      # Convert to matrix for algorithm
      X <- as.matrix(X_quanti)

      # -------------------------------------------------------------------------
      # AUTOMATIC K SELECTION (if n_clusters = "auto")
      # -------------------------------------------------------------------------
      if (is.character(private$.n_clusters) && private$.n_clusters == "auto") {

        message("Automatic K selection enabled. Testing K = ",
                min(private$.k_range), " to ", max(private$.k_range), "...")

        # Call the automatic K selection function from metrics_kmeans.R
        auto_result <- kmeans_find_optimal_k(
          data = X_quanti,
          k_range = private$.k_range,
          method = private$.selection_method,
          n_init = private$.n_init,
          random_state = private$.random_state,
          distance_metric = private$.distance_metric
        )

        # Extract optimal K
        optimal_k <- auto_result$optimal_k

        # Display informative message
        message(sprintf(
          "Optimal K selected: %d (method: %s, avg_silhouette: %.3f)",
          optimal_k,
          private$.selection_method,
          auto_result$metrics$avg_silhouette[auto_result$metrics$k == optimal_k]
        ))

        # Set n_clusters to the optimal value
        self$n_clusters <- optimal_k
      }

      # Check n_clusters doesn't exceed number of variables
      if (self$n_clusters > ncol(X)) {
        stop(paste0("n_clusters (", self$n_clusters,
                    ") cannot exceed number of variables (", ncol(X), ")"))
      }
      
      # Multiple initializations
      best_inertia <- -Inf
      best_result <- NULL

      for (init in 1:private$.n_init) {

        result <- private$kmeans_single_run(X)

        if (result$inertia > best_inertia) {
          best_inertia <- result$inertia
          best_result <- result
        }
      }
      
      # Store best results
      self$labels <- best_result$clusters
      private$.centroids <- best_result$centroids
      private$.inertia <- best_result$inertia
      private$.n_iter <- best_result$n_iter
      private$.cluster_inertias <- best_result$cluster_inertias
      
      # Mark as fitted
      self$fitted <- TRUE
      
      invisible(self)
    },
    
    # -----------------------------------------------------------------------
    # Predict method
    # -----------------------------------------------------------------------
    
    # Predict cluster assignments for new variables
    #
    # @description
    # Assigns new variables to existing clusters based on fitted centroids.
    #
    # @param data data.frame or matrix with same number of observations as training data.
    #
    # @return List with:
    # \itemize{
    #   \item clusters: Integer vector of cluster assignments
    #   \item distances: Matrix of distances to each centroid
    #   \item correlations: Matrix of correlations with each centroid
    # }
    predict = function(data) {
      
      if (!self$fitted) {
        stop("Model must be fitted before calling predict(). Use $fit() first.")
      }
      
      # Validate new data
      if (!is.data.frame(data) && !is.matrix(data)) {
        stop("'data' must be a data.frame or matrix")
      }
      
      X_new <- as.matrix(data)
      
      # Check number of observations matches
      if (nrow(X_new) != nrow(private$.centroids)) {
        stop(paste0("New data must have ", nrow(private$.centroids), 
                    " observations (same as training data)"))
      }
      
      # Compute distances and correlations
      result <- private$compute_distances(X_new, private$.centroids, return_correlations = TRUE)

      # Assign to nearest clusters
      clusters_new <- private$assign_clusters(result$distances)
      names(clusters_new) <- colnames(X_new)

      return(list(
        clusters = clusters_new,
        distances = result$distances,
        correlations = result$correlations
      ))
    },
    
    # -----------------------------------------------------------------------
    # Print method
    # -----------------------------------------------------------------------
    
    # Print method for KmeansVariables
    #
    # @description
    # Displays concise information about the fitted model.
    #
    print = function() {
      
      cat("K-means Variable Clustering\n")
      cat("===========================\n\n")
      
      if (!self$fitted) {
        cat("Model not fitted yet. Use $fit(data) to train the model.\n")
        return(invisible(self))
      }
      
      cat("Number of clusters:", self$n_clusters, "\n")
      cat("Number of variables:", length(self$labels), "\n")
      cat("Iterations:", private$.n_iter, "\n")
      cat("Total inertia:", round(private$.inertia, 4), "\n\n")
      
      cat("Cluster sizes:\n")
      print(table(self$labels))
      
      invisible(self)
    },
    
    # -----------------------------------------------------------------------
    # Summary method
    # -----------------------------------------------------------------------
    
    # Summary method for KmeansVariables
    #
    # @description
    # Displays detailed information about the fitted model.
    #
    summary = function() {
      
      if (!self$fitted) {
        cat("Model not fitted yet. Use $fit(data) to train the model.\n")
        return(invisible(self))
      }
      
      self$print()
      cat("\n")
      
      # Cluster-level statistics
      cat("Cluster Statistics:\n")
      cat("===================\n\n")
      
      X <- as.matrix(self$get_quanti_data())
      var_names <- colnames(X)
      
      for (k in 1:self$n_clusters) {
        
        var_indices <- which(self$labels == k)
        n_vars <- length(var_indices)
        
        cat(sprintf("Cluster %d (%d variables):\n", k, n_vars))
        cat(sprintf("  Within-cluster inertia: %.4f\n", 
                    private$.cluster_inertias[k]))
        
        if (n_vars > 0) {
          # Mean correlation with centroid
          cors <- sapply(var_indices, function(j) {
            cor(X[, j], private$.centroids[, k])
          })
          cat(sprintf("  Mean correlation with centroid: %.4f\n", 
                      mean(abs(cors))))
          
          # Top contributing variables
          contributions <- cors^2
          names(contributions) <- var_names[var_indices]
          top_vars <- head(sort(contributions, decreasing = TRUE), 3)
          
          cat("  Top variables:\n")
          for (i in seq_along(top_vars)) {
            cat(sprintf("    %s (R² = %.4f)\n", 
                        names(top_vars)[i], top_vars[i]))
          }
        }
        cat("\n")
      }

      invisible(self)
    },

    # -----------------------------------------------------------------------
    # METRIC WRAPPER METHODS
    # -----------------------------------------------------------------------

    # Calculate silhouette coefficient
    #
    # Wrapper for kmeans_silhouette() that uses the fitted model
    #
    # @return Data frame with silhouette values per variable
    silhouette = function() {
      if (!self$fitted) {
        stop("Model must be fitted before calculating silhouette. Use $fit() first.")
      }

      data <- self$get_quanti_data()
      return(kmeans_silhouette(data, self$labels, private$.centroids, private$.distance_metric))
    },

    # Calculate Calinski-Harabasz index
    #
    # Wrapper for kmeans_calinski_harabasz() that uses the fitted model
    #
    # @return Numeric Calinski-Harabasz index
    calinski_harabasz = function() {
      if (!self$fitted) {
        stop("Model must be fitted before calculating CH index. Use $fit() first.")
      }

      data <- self$get_quanti_data()
      return(kmeans_calinski_harabasz(data, self$labels, private$.centroids, private$.distance_metric))
    },

    # Calculate intra-cluster correlations
    #
    # Wrapper for kmeans_intra_correlation() that uses the fitted model
    #
    # @return Data frame with mean correlation per cluster
    intra_correlation = function() {
      if (!self$fitted) {
        stop("Model must be fitted before calculating intra-cluster correlation. Use $fit() first.")
      }

      data <- self$get_quanti_data()
      return(kmeans_intra_correlation(data, self$labels))
    },

    # Calculate variable contributions
    #
    # Wrapper for kmeans_contributions() that uses the fitted model
    #
    # @return Data frame with contributions per variable
    contributions = function() {
      if (!self$fitted) {
        stop("Model must be fitted before calculating contributions. Use $fit() first.")
      }

      data <- self$get_quanti_data()
      return(kmeans_contributions(data, self$labels, private$.centroids))
    },

    # Calculate correlation table
    #
    # Wrapper for kmeans_correlation_table() that uses the fitted model
    #
    # @param round_digits Number of decimal places to round to. Default: 3.
    # @return Data frame with variable-cluster correlations
    correlation_table = function(round_digits = 3) {
      if (!self$fitted) {
        stop("Model must be fitted before calculating correlation table. Use $fit() first.")
      }

      data <- self$get_quanti_data()
      return(kmeans_correlation_table(data, self$labels, private$.centroids, round_digits))
    },

    # -----------------------------------------------------------------------
    # PLOT WRAPPER METHODS
    # -----------------------------------------------------------------------

    # Plot elbow curve for K selection
    #
    # Wrapper for plot_kmeans_elbow()
    #
    # @param k_range Range of K values to test. Default: 2:10.
    # @param ... Additional arguments passed to plot_kmeans_elbow()
    # @return None (displays plot)
    plot_elbow = function(k_range = 2:10, ...) {
      data <- self$get_quanti_data()

      # Compute elbow data
      elbow_data <- kmeans_elbow(data, k_range = k_range,
                                 n_init = private$.n_init,
                                 random_state = private$.random_state,
                                 distance_metric = private$.distance_metric)

      # Plot with fitted K highlighted if model is fitted
      plot_kmeans_elbow(elbow_data, ...)

      # Add red point for current K if fitted
      if (self$fitted) {
        current_k <- self$n_clusters
        if (current_k %in% elbow_data$k) {
          points(current_k, elbow_data$inertia[elbow_data$k == current_k],
                 pch = 19, col = "red", cex = 2)
        }
      }
    },

    # Plot contributions
    #
    # Wrapper for plot_kmeans_contributions() that uses the fitted model
    #
    # @param top_n Number of top variables to display per cluster. Default: 5.
    # @param ... Additional arguments passed to plot_kmeans_contributions()
    # @return None (displays plot)
    plot_contributions = function(top_n = 5, ...) {
      if (!self$fitted) {
        stop("Model must be fitted before plotting. Use $fit() first.")
      }

      contrib_data <- self$contributions()
      plot_kmeans_contributions(contrib_data, top_n = top_n, ...)
    },

    # Plot variable projection (PCA-style)
    #
    # Wrapper for plot_kmeans_projection() that uses the fitted model
    #
    # @param supplementary_variables Optional data.frame/matrix of supplementary
    #   variables to project onto the PCA visualization
    # @param supplementary_clusters Optional cluster assignments for supplementary
    #   variables (from predict())
    # @param ... Additional arguments passed to plot_kmeans_projection()
    # @return None (displays plot)
    plot_projection = function(supplementary_variables = NULL,
                               supplementary_clusters = NULL,
                               ...) {
      if (!self$fitted) {
        stop("Model must be fitted before plotting. Use $fit() first.")
      }

      data <- self$get_quanti_data()
      plot_kmeans_projection(data, self$labels,
                            supplementary_variables = supplementary_variables,
                            supplementary_clusters = supplementary_clusters,
                            ...)
    }
  ),
  
  # ==========================================================================
  # ACTIVE BINDINGS (getters)
  # ==========================================================================
  
  active = list(
    
    # Returns cluster assignments for each variable
    clusters = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(self$labels)
    },
    
    # Returns matrix of cluster centroids (first PCs)
    centroids = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.centroids)
    },
    
    # Returns total within-cluster inertia
    inertia = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.inertia)
    },
    
    # Returns number of iterations until convergence
    n_iter = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.n_iter)
    },
    
    # Returns number of variables per cluster
    cluster_sizes = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      sizes <- table(self$labels)
      names(sizes) <- paste0("Cluster", 1:self$n_clusters)
      return(sizes)
    },
    
    # Returns within-cluster inertia for each cluster
    cluster_inertias = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.cluster_inertias)
    }
  )
)