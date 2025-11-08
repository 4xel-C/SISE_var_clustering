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
#' - Distance calculation based on squared correlation with centroids
#' - Centroid computation via Principal Component Analysis (PCA)
#' - Multiple random initializations to avoid local minima
#' - Convergence monitoring and iteration tracking
#'
#' @section Algorithm:
#' The K-means for variables works as follows:
#' \enumerate{
#'   \item **Initialization**: Randomly assign variables to K clusters
#'   \item **Centroid Calculation**: For each cluster, compute the first principal 
#'         component (via PCA) as the cluster centroid
#'   \item **Distance Calculation**: Compute distance between each variable and each 
#'         centroid using: \code{d(X_j, c_k) = sqrt(1 - cor(X_j, c_k)^2)}
#'   \item **Reallocation**: Reassign each variable to its nearest centroid
#'   \item **Convergence Check**: Repeat steps 2-4 until no reassignment occurs or 
#'         max iterations reached
#' }
#'
#' @section Distance Metric:
#' The distance between a variable \eqn{X_j} and a centroid \eqn{c_k} is defined as:
#' \deqn{d(X_j, c_k) = \sqrt{1 - \rho^2(X_j, c_k)}}
#' where \eqn{\rho(X_j, c_k)} is the Pearson correlation coefficient.
#' 
#' A high correlation (close to ±1) results in a small distance (close to 0),
#' indicating that the variable is well-represented by the centroid.
#'
#' @section Objective Function:
#' The algorithm minimizes the within-cluster inertia:
#' \deqn{W = \sum_{k=1}^{K} \sum_{X_j \in C_k} [1 - \rho^2(X_j, c_k)]}
#' where \eqn{C_k} denotes the set of variables in cluster \eqn{k}.
#'
#' @section Public Methods:
#' \describe{
#'   \item{\code{initialize(n_clusters = 3, max_iter = 100, tol = 1e-4, 
#'                          n_init = 10, random_state = NULL)}}{
#'     Creates a new instance of KmeansVariables.
#'     \itemize{
#'       \item \code{n_clusters}: Number of clusters (default: 3)
#'       \item \code{max_iter}: Maximum number of iterations (default: 100)
#'       \item \code{tol}: Convergence tolerance (default: 1e-4)
#'       \item \code{n_init}: Number of random initializations (default: 10)
#'       \item \code{random_state}: Seed for reproducibility (default: NULL)
#'     }
#'   }
#'   \item{\code{fit(data)}}{
#'     Fits the K-means clustering model on the provided data.
#'     \itemize{
#'       \item \code{data}: Data frame containing quantitative variables to cluster
#'     }
#'     This method:
#'     \itemize{
#'       \item Validates and loads the data (must have ≥2 quantitative variables)
#'       \item Runs the algorithm \code{n_init} times with different initializations
#'       \item Keeps the solution with the lowest inertia
#'       \item Stores cluster assignments, centroids, and performance metrics
#'     }
#'   }
#'   \item{\code{predict(data)}}{
#'     Assigns new variables to existing clusters.
#'     \itemize{
#'       \item \code{data}: Data frame with same number of observations as training data
#'     }
#'     Returns a list with:
#'     \itemize{
#'       \item \code{clusters}: Integer vector of cluster assignments
#'       \item \code{distances}: Matrix of distances to each centroid
#'       \item \code{correlations}: Matrix of correlations with each centroid
#'     }
#'   }
#'   \item{\code{print()}}{
#'     Displays concise information about the fitted model:
#'     number of clusters, iterations, and final inertia.
#'   }
#'   \item{\code{summary()}}{
#'     Displays detailed information including:
#'     \itemize{
#'       \item Cluster sizes (number of variables per cluster)
#'       \item Within-cluster inertia per cluster
#'       \item Mean correlation within each cluster
#'       \item Top contributing variables per cluster
#'     }
#'   }
#' }
#'
#' @section Active Bindings (getters):
#' \describe{
#'   \item{\code{clusters}}{
#'     Returns an integer vector of cluster assignments for each variable
#'   }
#'   \item{\code{centroids}}{
#'     Returns a matrix where each column is a cluster centroid 
#'     (first principal component)
#'   }
#'   \item{\code{inertia}}{
#'     Returns the total within-cluster inertia (lower is better)
#'   }
#'   \item{\code{n_iter}}{
#'     Returns the number of iterations performed until convergence
#'   }
#'   \item{\code{cluster_sizes}}{
#'     Returns a named vector with the number of variables in each cluster
#'   }
#'   \item{\code{cluster_inertias}}{
#'     Returns a named vector with the within-cluster inertia for each cluster
#'   }
#' }
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
#' To avoid local minima, the algorithm runs \code{n_init} times with different
#' random initializations. The solution with the lowest inertia is retained.
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
#' # Fit K-means with 3 clusters
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
#' @note
#' This class is not exported; it is intended for internal use
#' within the package.
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
#' @keywords internal
#'
#' @noRd
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
    
    # Results storage
    .centroids = NULL,       # Matrix of cluster centroids (PCA components)
    .inertia = NULL,         # Total within-cluster inertia
    .n_iter = NULL,          # Number of iterations performed
    .cluster_inertias = NULL, # Inertia per cluster
    
    # ==========================================================================
    # PRIVATE METHODS - CORE ALGORITHM
    # ==========================================================================
    
    #' Validate algorithm parameters
    #' 
    #' Checks that max_iter, tol, n_init are valid positive values
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
    
    #' Calculate centroid of a cluster using PCA
    #' 
    #' Computes the first principal component of variables in a cluster.
    #' This component serves as the cluster centroid.
    #' 
    #' @param X Matrix or data.frame of observations × variables
    #' @param var_indices Integer vector of column indices in the cluster
    #' @return Numeric vector (centroid) of length nrow(X)
    calculate_centroid = function(X, var_indices) {
      
      # Handle single variable case
      if (length(var_indices) == 1) {
        return(as.numeric(scale(X[, var_indices])))
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
          return(as.numeric(scale(X_cluster[, 1])))
        } else {
          # No variables left, return zeros
          return(rep(0, nrow(X_cluster)))
        }
      }
      
      # Perform PCA (scale variables first)
      pca_result <- tryCatch({
        prcomp(X_cluster, center = TRUE, scale. = TRUE)
      }, error = function(e) {
        # If PCA fails, return scaled mean
        return(list(x = cbind(scale(rowMeans(X_cluster)))))
      })
      
      # Return first principal component (scores)
      centroid <- pca_result$x[, 1]
      
      # Ensure positive direction for consistency
      if (sum(centroid) < 0) {
        centroid <- -centroid
      }
      
      return(as.numeric(centroid))
    },
    
    #' Compute distances between variables and centroids
    #' 
    #' Calculates the distance matrix based on squared correlation.
    #' Distance = sqrt(1 - cor(X_j, c_k)^2)
    #' 
    #' @param X Matrix or data.frame (n observations × p variables)
    #' @param centroids Matrix (n observations × K clusters)
    #' @return Matrix (p variables × K clusters) of distances
    compute_distances = function(X, centroids) {

      p <- ncol(X)
      K <- ncol(centroids)

      dist_matrix <- matrix(0, nrow = p, ncol = K)
      rownames(dist_matrix) <- colnames(X)
      colnames(dist_matrix) <- paste0("Cluster", 1:K)

      for (j in 1:p) {
        for (k in 1:K) {
          # Correlation between variable j and centroid k
          cor_val <- suppressWarnings(cor(X[, j], centroids[, k]))

          # Handle NA from zero variance: set distance to 1 (maximum)
          if (is.na(cor_val)) {
            cor_val <- 0
          }

          # Distance = sqrt(1 - cor^2)
          dist_matrix[j, k] <- sqrt(1 - cor_val^2)
        }
      }

      return(dist_matrix)
    },
    
    #' Assign variables to nearest centroids
    #' 
    #' @param dist_matrix Matrix (p variables × K clusters) of distances
    #' @return Integer vector of cluster assignments (1 to K)
    assign_clusters = function(dist_matrix) {
      
      # For each variable, find cluster with minimum distance
      clusters <- apply(dist_matrix, 1, which.min)
      
      return(as.integer(clusters))
    },
    
    #' Calculate total within-cluster inertia
    #' 
    #' Inertia = sum over all variables of (1 - cor^2 with assigned centroid)
    #' 
    #' @param X Matrix or data.frame (n × p)
    #' @param clusters Integer vector of cluster assignments
    #' @param centroids Matrix (n × K) of centroids
    #' @return Numeric value (total inertia)
    calculate_inertia = function(X, clusters, centroids) {

      p <- ncol(X)
      inertia <- 0

      for (j in 1:p) {
        k <- clusters[j]  # Assigned cluster
        cor_val <- suppressWarnings(cor(X[, j], centroids[, k]))

        # Handle NA from zero variance: set correlation to 0
        if (is.na(cor_val)) {
          cor_val <- 0
        }

        inertia <- inertia + (1 - cor_val^2)
      }

      return(inertia)
    },
    
    #' Calculate within-cluster inertia for each cluster
    #' 
    #' @param X Matrix or data.frame (n × p)
    #' @param clusters Integer vector of cluster assignments
    #' @param centroids Matrix (n × K) of centroids
    #' @return Named numeric vector (inertia per cluster)
    calculate_cluster_inertias = function(X, clusters, centroids) {
      
      K <- ncol(centroids)
      cluster_inertias <- numeric(K)
      names(cluster_inertias) <- paste0("Cluster", 1:K)
      
      for (k in 1:K) {
        var_indices <- which(clusters == k)
        
        if (length(var_indices) == 0) {
          cluster_inertias[k] <- 0
          next
        }
        
        inertia_k <- 0
        for (j in var_indices) {
          cor_val <- suppressWarnings(cor(X[, j], centroids[, k]))

          # Handle NA from zero variance: set correlation to 0
          if (is.na(cor_val)) {
            cor_val <- 0
          }

          inertia_k <- inertia_k + (1 - cor_val^2)
        }

        cluster_inertias[k] <- inertia_k
      }
      
      return(cluster_inertias)
    },
    
    #' Single run of K-means algorithm
    #' 
    #' @param X Matrix or data.frame (n × p)
    #' @return List with clusters, centroids, inertia, n_iter
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
      
      inertia_old <- Inf
      
      for (iter in 1:private$.max_iter) {
        
        # Step 1: Calculate centroids
        centroids <- matrix(0, nrow = n, ncol = K)
        colnames(centroids) <- paste0("Cluster", 1:K)
        
        for (k in 1:K) {
          var_indices <- which(clusters == k)
          centroids[, k] <- private$calculate_centroid(X, var_indices)
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
        
        # Step 4: Calculate inertia
        inertia_new <- private$calculate_inertia(X, clusters_new, centroids)

        # Check convergence
        # Use isTRUE(all.equal()) to handle NA values
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
      
      # Calculate final metrics
      final_inertia <- private$calculate_inertia(X, clusters, centroids)
      cluster_inertias <- private$calculate_cluster_inertias(X, clusters, centroids)
      
      return(list(
        clusters = clusters,
        centroids = centroids,
        inertia = final_inertia,
        n_iter = iter,
        cluster_inertias = cluster_inertias
      ))
    }
  ),
  
  # ==========================================================================
  # PUBLIC METHODS
  # ==========================================================================
  
  public = list(
    
    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------
    
    #' @title Constructor for KmeansVariables
    #'
    #' @description
    #' Initializes a new instance of the `KmeansVariables` class.
    #'
    #' @param n_clusters Integer. Number of clusters. Default: 3.
    #' @param max_iter Integer. Maximum number of iterations. Default: 100.
    #' @param tol Numeric. Convergence tolerance. Default: 1e-4.
    #' @param n_init Integer. Number of random initializations. Default: 10.
    #' @param random_state Integer or NULL. Random seed for reproducibility. Default: NULL.
    #'
    #' @return A new instance of `KmeansVariables`.
    #' @noRd
    initialize = function(n_clusters = 3, max_iter = 100, tol = 1e-4, 
                          n_init = 10, random_state = NULL) {
      
      # Validate parameters
      private$validate_params(max_iter, tol, n_init)
      
      # Set n_clusters via parent class setter (includes validation)
      self$n_clusters <- n_clusters
      
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
    
    #' Fit K-means clustering on variables
    #'
    #' @description
    #' Performs K-means clustering on quantitative variables.
    #' Runs multiple initializations and keeps the best solution.
    #'
    #' @param data data.frame or matrix containing quantitative variables to cluster.
    #'
    #' @details
    #' - Validates that data contains at least 2 quantitative variables
    #' - If the dataset contains qualitative variables, they are automatically ignored
    #'   (a message is displayed to inform the user)
    #' - Runs the algorithm `n_init` times with different random starts
    #' - Keeps the solution with the lowest inertia
    #' - Stores results in private fields accessible via active bindings
    #'
    #' @return None. Results are stored internally and accessible via getters.
    #' @noRd
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

      # Check n_clusters doesn't exceed number of variables
      if (self$n_clusters > ncol(X)) {
        stop(paste0("n_clusters (", self$n_clusters, 
                    ") cannot exceed number of variables (", ncol(X), ")"))
      }
      
      # Multiple initializations
      best_inertia <- Inf
      best_result <- NULL
      
      for (init in 1:private$.n_init) {
        
        result <- private$kmeans_single_run(X)
        
        if (result$inertia < best_inertia) {
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
    
    #' Predict cluster assignments for new variables
    #'
    #' @description
    #' Assigns new variables to existing clusters based on fitted centroids.
    #'
    #' @param data data.frame or matrix with same number of observations as training data.
    #'
    #' @return List with:
    #' \itemize{
    #'   \item clusters: Integer vector of cluster assignments
    #'   \item distances: Matrix of distances to each centroid
    #'   \item correlations: Matrix of correlations with each centroid
    #' }
    #' @noRd
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
      
      # Compute distances to centroids
      dist_matrix <- private$compute_distances(X_new, private$.centroids)
      
      # Assign to nearest clusters
      clusters_new <- private$assign_clusters(dist_matrix)
      names(clusters_new) <- colnames(X_new)
      
      # Calculate correlations for interpretation
      K <- ncol(private$.centroids)
      cor_matrix <- matrix(0, nrow = ncol(X_new), ncol = K)
      rownames(cor_matrix) <- colnames(X_new)
      colnames(cor_matrix) <- paste0("Cluster", 1:K)
      
      for (j in 1:ncol(X_new)) {
        for (k in 1:K) {
          cor_matrix[j, k] <- cor(X_new[, j], private$.centroids[, k])
        }
      }
      
      return(list(
        clusters = clusters_new,
        distances = dist_matrix,
        correlations = cor_matrix
      ))
    },
    
    # -----------------------------------------------------------------------
    # Print method
    # -----------------------------------------------------------------------
    
    #' Print method for KmeansVariables
    #'
    #' @description
    #' Displays concise information about the fitted model.
    #'
    #' @noRd
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
    
    #' Summary method for KmeansVariables
    #'
    #' @description
    #' Displays detailed information about the fitted model.
    #'
    #' @noRd
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
    }
  ),
  
  # ==========================================================================
  # ACTIVE BINDINGS (getters)
  # ==========================================================================
  
  active = list(
    
    #' Returns cluster assignments for each variable
    clusters = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(self$labels)
    },
    
    #' Returns matrix of cluster centroids (first PCs)
    centroids = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.centroids)
    },
    
    #' Returns total within-cluster inertia
    inertia = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.inertia)
    },
    
    #' Returns number of iterations until convergence
    n_iter = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.n_iter)
    },
    
    #' Returns number of variables per cluster
    cluster_sizes = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      sizes <- table(self$labels)
      names(sizes) <- paste0("Cluster", 1:self$n_clusters)
      return(sizes)
    },
    
    #' Returns within-cluster inertia for each cluster
    cluster_inertias = function() {
      if (!self$fitted) {
        warning("Model not fitted. Returning NULL.")
        return(NULL)
      }
      return(private$.cluster_inertias)
    }
  )
)