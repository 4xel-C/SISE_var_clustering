# ===========================================================================
# Metrics for K-means Variable Clustering
# ===========================================================================
#
# Collection of functions to evaluate K-means variable clustering results:
# - kmeans_elbow: Elbow method (within-cluster homogeneity)
# - kmeans_silhouette: Silhouette coefficient
# - kmeans_calinski_harabasz: Calinski-Harabasz index
# - kmeans_intra_correlation: Within-cluster correlation analysis
# - kmeans_contributions: Variable contribution (R²) to clusters
# - kmeans_correlation_table: Correlation table variable-cluster
# - kmeans_find_optimal_k: Find optimal K using multiple criteria
#
# All functions are internal utilities (@noRd).
# ===========================================================================


#' Calculate Elbow Method Curve
#'
#' @description
#' Computes the within-cluster inertia for different values of K.
#' Used to identify the "elbow point" where adding more clusters
#' provides diminishing returns.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param k_range Integer vector of K values to test (e.g., 2:10)
#' @param n_init Integer. Number of random initializations per K. Default: 10.
#' @param random_state Integer or NULL. Random seed for reproducibility.
#' @param max_iter Integer. Maximum iterations for K-means. Default: 100.
#' @param tol Numeric. Convergence tolerance. Default: 1e-4.
#' @param distance_metric Character. Distance metric: "r_squared" or "r_signed". Default: "r_squared".
#'
#' @return data.frame with columns:
#' \itemize{
#'   \item k: Number of clusters
#'   \item inertia: Within-cluster inertia
#'   \item inertia_pct: Percentage of total variance explained
#' }
#'
#' @examples
#' \dontrun{
#' elbow_data <- kmeans_elbow(data, k_range = 2:8)
#' plot(elbow_data$k, elbow_data$inertia, type = "b")
#' }
#'
#' @noRd
kmeans_elbow <- function(data, k_range = 2:8, n_init = 10,
                         random_state = NULL, max_iter = 50, tol = 1e-3,
                         distance_metric = "r_squared") {

  if (!is.numeric(k_range) || any(k_range < 1)) {
    stop("'k_range' must be a vector of positive integers")
  }

  # Calculate total variance (for k=1, all variables in one cluster)
  X <- as.matrix(data)
  total_var <- sum(apply(as.matrix(X), 2, function(x) var(x)))

  results <- data.frame(
    k = integer(),
    inertia = numeric(),
    inertia_pct = numeric()
  )

  for (k in k_range) {

    # Fit K-means with k clusters and specified distance_metric
    km <- KmeansVariables$new(
      n_clusters = k,
      max_iter = max_iter,
      tol = tol,
      n_init = n_init,
      random_state = random_state,
      distance_metric = distance_metric
    )

    km$fit(data)

    # Store results
    results <- rbind(results, data.frame(
      k = k,
      inertia = km$inertia,
      inertia_pct = (1 - km$inertia / ncol(X)) * 100
    ))
  }

  return(results)
}


#' Calculate Silhouette Coefficient for Variable Clustering
#'
#' @description
#' Computes the silhouette coefficient for each variable, measuring
#' how well each variable fits within its assigned cluster.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param clusters Integer vector of cluster assignments
#' @param centroids Matrix of cluster centroids (n_obs × n_clusters)
#' @param distance_metric Character. Distance metric: "r_squared" or "r_signed". Default: "r_squared".
#'
#' @return data.frame with columns:
#' \itemize{
#'   \item variable: Variable name
#'   \item cluster: Assigned cluster
#'   \item silhouette: Silhouette coefficient (-1 to 1)
#' }
#'
#' @details
#' Silhouette coefficient for variable j:
#' \deqn{s(j) = \frac{b(j) - a(j)}{\max(a(j), b(j))}}
#' where:
#' \itemize{
#'   \item a(j) = average distance to variables in same cluster
#'   \item b(j) = average distance to variables in nearest other cluster
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item s(j) close to 1: well-clustered
#'   \item s(j) close to 0: on cluster boundary
#'   \item s(j) negative: possibly misclassified
#' }
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' sil <- kmeans_silhouette(data, km$clusters, km$centroids)
#' mean(sil$silhouette)  # Average silhouette score
#' }
#'
#' @noRd
kmeans_silhouette <- function(data, clusters, centroids, distance_metric = "r_squared") {

  X <- as.matrix(data)
  p <- ncol(X)
  K <- ncol(centroids)

  # Validate distance_metric
  if (!distance_metric %in% c("r_squared", "r_signed")) {
    stop("'distance_metric' must be either 'r_squared' or 'r_signed'")
  }

  # Compute distance matrix between all variables and all centroids
  dist_matrix <- matrix(0, nrow = p, ncol = K)

  for (j in 1:p) {
    for (k in 1:K) {
      cor_val <- cor(X[, j], centroids[, k])

      # Calculate distance based on distance_metric
      if (distance_metric == "r_squared") {
        dist_matrix[j, k] <- sqrt(1 - cor_val^2)
      } else {
        dist_matrix[j, k] <- sqrt(1 - cor_val)
      }
    }
  }

  silhouettes <- numeric(p)

  # Early check for single cluster case
  if (K == 1) {
    return(data.frame(
      variable = colnames(X),
      cluster = clusters,
      silhouette = rep(0, p)
    ))
  }

  for (j in 1:p) {

    k_own <- clusters[j]

    # a(j): average distance to own cluster centroid
    a_j <- dist_matrix[j, k_own]

    # b(j): distance to nearest other cluster centroid
    other_clusters <- setdiff(1:K, k_own)
    b_j <- min(dist_matrix[j, other_clusters])

    # Silhouette coefficient
    if (max(a_j, b_j) == 0) {
      silhouettes[j] <- 0
    } else {
      silhouettes[j] <- (b_j - a_j) / max(a_j, b_j)
    }
  }

  results <- data.frame(
    variable = colnames(X),
    cluster = clusters,
    silhouette = silhouettes
  )

  return(results)
}


#' Calculate Average Silhouette for Different K Values
#'
#' @description
#' Computes the average silhouette coefficient for a range of K values.
#' Helps identify the optimal number of clusters.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param k_range Integer vector of K values to test
#' @param n_init Integer. Number of random initializations per K. Default: 10.
#' @param random_state Integer or NULL. Random seed for reproducibility.
#' @param max_iter Integer. Maximum iterations. Default: 100.
#' @param tol Numeric. Convergence tolerance. Default: 1e-4.
#' @param distance_metric Character. Distance metric: "r_squared" or "r_signed". Default: "r_squared".
#'
#' @return data.frame with columns:
#' \itemize{
#'   \item k: Number of clusters
#'   \item avg_silhouette: Average silhouette coefficient
#'   \item min_silhouette: Minimum silhouette (worst variable)
#'   \item max_silhouette: Maximum silhouette (best variable)
#' }
#'
#' @examples
#' \dontrun{
#' sil_data <- kmeans_silhouette_range(data, k_range = 2:8)
#' # Optimal K = argmax(avg_silhouette)
#' best_k <- sil_data$k[which.max(sil_data$avg_silhouette)]
#' }
#'
#' @noRd
kmeans_silhouette_range <- function(data, k_range = 2:10, n_init = 10,
                                    random_state = NULL, max_iter = 100,
                                    tol = 1e-4, distance_metric = "r_squared") {

  results <- data.frame(
    k = integer(),
    avg_silhouette = numeric(),
    min_silhouette = numeric(),
    max_silhouette = numeric()
  )

  for (k in k_range) {

    # Fit K-means with specified distance_metric
    km <- KmeansVariables$new(
      n_clusters = k,
      max_iter = max_iter,
      tol = tol,
      n_init = n_init,
      random_state = random_state,
      distance_metric = distance_metric
    )

    km$fit(data)

    # Calculate silhouette with same distance_metric
    sil <- kmeans_silhouette(data, km$clusters, km$centroids, distance_metric)

    results <- rbind(results, data.frame(
      k = k,
      avg_silhouette = mean(sil$silhouette),
      min_silhouette = min(sil$silhouette),
      max_silhouette = max(sil$silhouette)
    ))
  }

  return(results)
}


#' Calculate Calinski-Harabasz Index
#'
#' @description
#' Computes the Calinski-Harabasz (CH) index, which measures
#' the ratio of between-cluster to within-cluster variance.
#' Higher values indicate better-defined clusters.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param clusters Integer vector of cluster assignments
#' @param centroids Matrix of cluster centroids
#' @param distance_metric Character. Distance metric used for clustering (for consistency checking).
#'
#' @return Numeric value (CH index). Higher is better.
#'
#' @details
#' The Calinski-Harabasz index is defined as:
#' \deqn{CH(K) = \frac{B/(K-1)}{W/(p-K)}}
#' where:
#' \itemize{
#'   \item B = between-cluster variance
#'   \item W = within-cluster variance (inertia)
#'   \item K = number of clusters
#'   \item p = number of variables
#' }
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' ch_score <- kmeans_calinski_harabasz(data, km$clusters, km$centroids)
#' }
#'
#' @noRd
kmeans_calinski_harabasz <- function(data, clusters, centroids, distance_metric = "r_squared") {

  X <- as.matrix(data)
  p <- ncol(X)
  K <- ncol(centroids)

  if (K == 1 || K == p) {
    warning("CH index undefined for K=1 or K=p")
    return(NA)
  }

  # Validate distance_metric
  if (!distance_metric %in% c("r_squared", "r_signed")) {
    stop("'distance_metric' must be either 'r_squared' or 'r_signed'")
  }

  # Within-cluster variance (W)
  # W = p - sum(eigenvalues) = p - total inertia
  W <- 0
  for (k in 1:K) {
    var_indices <- which(clusters == k)
    if (length(var_indices) == 0) next

    X_cluster <- X[, var_indices, drop = FALSE]

    # Single variable cluster
    if (length(var_indices) == 1) {
      W <- W + 0  # Perfect fit, variance = 0
      next
    }

    # Run PCA to get eigenvalue
    pca_result <- prcomp(X_cluster, center = TRUE, scale. = TRUE)
    eigenvalue <- pca_result$sdev[1]^2

    # Within-cluster variance = number of variables - eigenvalue
    W <- W + (length(var_indices) - eigenvalue)
  }

  # Between-cluster variance (B)
  B <- p - W

  # Calinski-Harabasz index
  CH <- (B / (K - 1)) / (W / (p - K))

  return(CH)
}


#' Calculate Calinski-Harabasz Index for Range of K
#'
#' @description
#' Computes CH index for different values of K to help identify
#' the optimal number of clusters.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param k_range Integer vector of K values to test
#' @param n_init Integer. Number of random initializations. Default: 10.
#' @param random_state Integer or NULL. Random seed.
#' @param max_iter Integer. Maximum iterations. Default: 100.
#' @param tol Numeric. Convergence tolerance. Default: 1e-4.
#' @param distance_metric Character. Distance metric: "r_squared" or "r_signed". Default: "r_squared".
#'
#' @return data.frame with columns:
#' \itemize{
#'   \item k: Number of clusters
#'   \item ch_index: Calinski-Harabasz index
#' }
#'
#' @examples
#' \dontrun{
#' ch_data <- kmeans_calinski_harabasz_range(data, k_range = 2:8)
#' # Optimal K = argmax(ch_index)
#' best_k <- ch_data$k[which.max(ch_data$ch_index)]
#' }
#'
#' @noRd
kmeans_calinski_harabasz_range <- function(data, k_range = 2:10, n_init = 10,
                                           random_state = NULL, max_iter = 100,
                                           tol = 1e-4, distance_metric = "r_squared") {

  results <- data.frame(
    k = integer(),
    ch_index = numeric()
  )

  for (k in k_range) {

    # Skip K=1 (undefined)
    if (k == 1) next

    # Fit K-means with specified distance_metric
    km <- KmeansVariables$new(
      n_clusters = k,
      max_iter = max_iter,
      tol = tol,
      n_init = n_init,
      random_state = random_state,
      distance_metric = distance_metric
    )

    km$fit(data)

    # Calculate CH index with same distance_metric
    ch <- kmeans_calinski_harabasz(data, km$clusters, km$centroids, distance_metric)

    results <- rbind(results, data.frame(
      k = k,
      ch_index = ch
    ))
  }

  return(results)
}


#' Calculate Within-Cluster Correlation Matrix
#'
#' @description
#' Computes the average correlation between variables within each cluster.
#' High intra-cluster correlation indicates cohesive clusters.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param clusters Integer vector of cluster assignments
#'
#' @return Named numeric vector with average correlation per cluster
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' intra_cor <- kmeans_intra_correlation(data, km$clusters)
#' print(intra_cor)
#' }
#'
#' @noRd
kmeans_intra_correlation <- function(data, clusters) {

  X <- as.matrix(data)
  K <- max(clusters)

  intra_cors <- numeric(K)
  names(intra_cors) <- paste0("Cluster", 1:K)

  for (k in 1:K) {

    var_indices <- which(clusters == k)
    n_vars <- length(var_indices)

    if (n_vars < 2) {
      intra_cors[k] <- NA
      next
    }

    # Use eigenvalue to measure cohesion
    X_k <- X[, var_indices, drop = FALSE]
    pca_result <- prcomp(X_k, center = TRUE, scale. = TRUE)
    eigenvalue <- pca_result$sdev[1]^2

    # Cohesion = eigenvalue / n_vars (proportion of variance explained by PC1)
    intra_cors[k] <- eigenvalue / n_vars
  }

  return(intra_cors)
}


#' Calculate Variable Contributions to Clusters
#'
#' @description
#' Computes the contribution (squared correlation) of each variable
#' to its assigned cluster centroid.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param clusters Integer vector of cluster assignments
#' @param centroids Matrix of cluster centroids
#'
#' @return data.frame with columns:
#' \itemize{
#'   \item variable: Variable name
#'   \item cluster: Assigned cluster
#'   \item contribution: R² with centroid (0 to 1)
#'   \item correlation: Pearson correlation with centroid
#' }
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' contrib <- kmeans_contributions(data, km$clusters, km$centroids)
#'
#' # Top 5 contributors per cluster
#' top_contrib <- contrib %>%
#'   group_by(cluster) %>%
#'   arrange(desc(contribution)) %>%
#'   slice_head(n = 5)
#' }
#'
#' @noRd
kmeans_contributions <- function(data, clusters, centroids) {

  X <- as.matrix(data)
  p <- ncol(X)

  contributions <- data.frame(
    variable = colnames(X),
    cluster = clusters,
    contribution = numeric(p),
    correlation = numeric(p)
  )

  for (j in 1:p) {
    k <- clusters[j]
    cor_val <- cor(X[, j], centroids[, k])
    contributions$correlation[j] <- cor_val
    # Contribution = R² = squared correlation = cos² (angle with PC1)
    contributions$contribution[j] <- cor_val^2
  }

  return(contributions)
}


#' Calculate Between/Within Cluster Ratio
#'
#' @description
#' Computes the ratio of between-cluster variance to within-cluster variance.
#' Higher values indicate better cluster separation.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param clusters Integer vector of cluster assignments
#' @param centroids Matrix of cluster centroids
#' @param distance_metric Character. Distance metric used for clustering (for consistency checking).
#'
#' @return Numeric value (ratio B/W). Higher is better.
#'
#' @details
#' This ratio measures cluster quality:
#' \itemize{
#'   \item High ratio: clusters are well-separated
#'   \item Low ratio: clusters overlap
#' }
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' ratio <- kmeans_between_within_ratio(data, km$clusters, km$centroids)
#' }
#'
#' @noRd
kmeans_between_within_ratio <- function(data, clusters, centroids, distance_metric = "r_squared") {

  X <- as.matrix(data)
  p <- ncol(X)
  K <- max(clusters)

  # Validate distance_metric
  if (!distance_metric %in% c("r_squared", "r_signed")) {
    stop("'distance_metric' must be either 'r_squared' or 'r_signed'")
  }

  # Within-cluster variance using eigenvalues
  W <- 0
  for (k in 1:K) {
    var_indices <- which(clusters == k)
    if (length(var_indices) == 0) next

    X_cluster <- X[, var_indices, drop = FALSE]

    # Single variable cluster
    if (length(var_indices) == 1) {
      W <- W + 0
      next
    }

    # Run PCA to get eigenvalue
    pca_result <- prcomp(X_cluster, center = TRUE, scale. = TRUE)
    eigenvalue <- pca_result$sdev[1]^2

    # Within-cluster variance = number of variables - eigenvalue
    W <- W + (length(var_indices) - eigenvalue)
  }

  # Between-cluster variance
  B <- p - W

  # Ratio
  if (W == 0) {
    return(Inf)
  }

  return(B / W)
}


#' Find Knee/Elbow Point in a Curve
#'
#' @description
#' Finds the elbow point in a curve using the maximum perpendicular distance method.
#' Draws a line from first to last point and finds the point with maximum distance to this line.
#'
#' @param x Numeric vector of x values (e.g., k values)
#' @param y Numeric vector of y values (e.g., metric values)
#' @param direction Character. Direction of the elbow: "decreasing" for metrics that decrease
#'        with K (like inertia), "increasing" for metrics that increase (like silhouette).
#'        Default: "decreasing".
#'
#' @return Integer index of the knee point in the input vectors
#'
#' @noRd
find_knee_point <- function(x, y, direction = "decreasing") {
  n <- length(x)
  if (n < 3) return(1)

  # Normalize x and y to [0, 1] for fair distance calculation
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))

  # Line from first to last point: ax + by + c = 0
  # Points: (x_norm[1], y_norm[1]) and (x_norm[n], y_norm[n])
  a <- y_norm[n] - y_norm[1]
  b <- x_norm[1] - x_norm[n]
  c <- x_norm[n] * y_norm[1] - x_norm[1] * y_norm[n]

  # Calculate perpendicular distance from each point to the line
  distances <- abs(a * x_norm + b * y_norm + c) / sqrt(a^2 + b^2)

  # For increasing curves, we want the point above the line with max distance
  # For decreasing curves, we want the point below the line with max distance
  if (direction == "increasing") {
    # For increasing metrics, the elbow is where curve is above the line
    above_line <- (a * x_norm + b * y_norm + c) > 0
    distances[!above_line] <- 0
  } else {
    # For decreasing metrics, the elbow is where curve is below the line
    below_line <- (a * x_norm + b * y_norm + c) < 0
    distances[!below_line] <- 0
  }

  knee_idx <- which.max(distances)
  return(knee_idx)
}


#' Find Optimal Number of Clusters
#'
#' @description
#' Automatically determines the optimal K using the elbow/knee detection method.
#' Uses perpendicular distance to find where marginal gains diminish.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param k_range Integer vector of K values to test. Default: 2:10.
#' @param method Character. Method to use: "silhouette", "calinski", or "all".
#'        Default: "all".
#' @param n_init Integer. Number of initializations. Default: 10.
#' @param random_state Integer or NULL. Random seed.
#' @param distance_metric Character. Distance metric: "r_squared" or "r_signed". Default: "r_squared".
#'
#' @return List with:
#' \itemize{
#'   \item optimal_k: Recommended number of clusters
#'   \item metrics: data.frame with all computed metrics
#'   \item method: Method used for recommendation
#' }
#'
#' @examples
#' \dontrun{
#' result <- kmeans_find_optimal_k(data, k_range = 2:8)
#' print(result$optimal_k)
#' print(result$metrics)
#' }
#'
#' @noRd
kmeans_find_optimal_k <- function(data, k_range = 2:10, method = "all",
                                  n_init = 10, random_state = NULL,
                                  distance_metric = "r_squared") {

  if (!method %in% c("silhouette", "calinski", "all")) {
    stop("'method' must be 'silhouette', 'calinski', or 'all'")
  }

  # Compute metrics with specified distance_metric
  elbow_data <- kmeans_elbow(data, k_range, n_init, random_state,
                            distance_metric = distance_metric)
  sil_data <- kmeans_silhouette_range(data, k_range, n_init, random_state,
                                      distance_metric = distance_metric)
  ch_data <- kmeans_calinski_harabasz_range(data, k_range, n_init, random_state,
                                            distance_metric = distance_metric)

  # Merge results
  metrics <- merge(elbow_data, sil_data, by = "k")
  metrics <- merge(metrics, ch_data, by = "k")

  # Determine optimal K using knee point detection
  if (method == "silhouette") {
    # Silhouette increases then decreases - find the knee where gains diminish
    knee_idx <- find_knee_point(metrics$k, metrics$avg_silhouette, direction = "increasing")
    optimal_k <- metrics$k[knee_idx]
  } else if (method == "calinski") {
    # Calinski-Harabasz increases then decreases - find the knee
    knee_idx <- find_knee_point(metrics$k, metrics$ch_index, direction = "increasing")
    optimal_k <- metrics$k[knee_idx]
  } else {
    # Combine silhouette and calinski methods
    knee_sil <- find_knee_point(metrics$k, metrics$avg_silhouette, direction = "increasing")
    knee_ch <- find_knee_point(metrics$k, metrics$ch_index, direction = "increasing")

    k_sil <- metrics$k[knee_sil]
    k_ch <- metrics$k[knee_ch]

    if (k_sil == k_ch) {
      optimal_k <- k_sil
    } else {
      # If disagreement, prefer silhouette
      optimal_k <- k_sil
      message("Methods disagree (silhouette: ", k_sil, ", calinski: ", k_ch,
              "). Using silhouette criterion.")
    }
  }

  return(list(
    optimal_k = optimal_k,
    metrics = metrics,
    method = method
  ))
}


#' Create Correlation Table: Variables vs Clusters
#'
#' @description
#' Creates a comprehensive table showing the correlation of each variable
#' with each cluster centroid. Helps identify the strength of cluster membership
#' and detect variables that might be misclassified or ambiguous.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param clusters Integer vector of cluster assignments
#' @param centroids Matrix of cluster centroids
#' @param round_digits Integer. Number of decimals to round correlations. Default: 3
#' @param distance_metric Character. Distance metric used for clustering: "r_squared" or "r_signed".
#'        Determines whether to use absolute correlation (r_squared) or signed correlation (r_signed)
#'        for calculating max_cor and separation metrics. Default: "r_squared".
#'
#' @return data.frame with columns:
#' \itemize{
#'   \item variable: Variable name
#'   \item assigned_cluster: Cluster to which the variable is assigned
#'   \item Cluster1, Cluster2, ...: Correlation with each cluster centroid
#'   \item max_cor: Maximum absolute correlation (with assigned cluster)
#'   \item second_max_cor: Second highest absolute correlation
#'   \item separation: Difference between max and second_max (higher = better separation)
#' }
#'
#' @details
#' This table is particularly useful to:
#' \itemize{
#'   \item Verify cluster assignments (high correlation with assigned cluster)
#'   \item Detect ambiguous variables (similar correlations with multiple clusters)
#'   \item Identify outlier variables (low correlation with all clusters)
#'   \item Assess cluster quality (high separation values indicate clear assignments)
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item \code{max_cor > 0.7}: Strong cluster membership
#'   \item \code{separation > 0.3}: Clear cluster assignment
#'   \item \code{separation < 0.1}: Ambiguous assignment (variable between clusters)
#' }
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' cor_table <- kmeans_correlation_table(data, km$clusters, km$centroids)
#'
#' # View full table
#' print(cor_table)
#'
#' # Find ambiguous variables (low separation)
#' ambiguous <- cor_table[cor_table$separation < 0.1, ]
#'
#' # Find strongly assigned variables
#' strong <- cor_table[cor_table$max_cor > 0.7, ]
#' }
#'
#' @noRd
kmeans_correlation_table <- function(data, clusters, centroids, round_digits = 3,
                                     distance_metric = "r_squared") {

  X <- as.matrix(data)
  p <- ncol(X)
  K <- ncol(centroids)
  var_names <- colnames(X)

  # Validate distance_metric
  if (!distance_metric %in% c("r_squared", "r_signed")) {
    stop("'distance_metric' must be either 'r_squared' or 'r_signed'")
  }

  # Initialize result data frame
  result <- data.frame(
    variable = var_names,
    assigned_cluster = clusters,
    stringsAsFactors = FALSE
  )

  # Calculate correlation of each variable with each cluster centroid
  cor_matrix <- matrix(NA, nrow = p, ncol = K)
  colnames(cor_matrix) <- paste0("Cluster", 1:K)

  for (j in 1:p) {
    for (k in 1:K) {
      cor_val <- suppressWarnings(cor(X[, j], centroids[, k]))

      # Handle NA from zero variance
      if (is.na(cor_val)) {
        cor_val <- 0
      }

      cor_matrix[j, k] <- cor_val
    }
  }

  # Add correlation columns to result
  result <- cbind(result, round(cor_matrix, round_digits))

  # Calculate additional metrics based on distance_metric
  # For r_squared: use absolute correlation (considers both positive and negative as similar)
  # For r_signed: use raw signed correlation (only positive values indicate similarity)
  if (distance_metric == "r_squared") {
    metric_matrix <- abs(cor_matrix)
  } else {
    metric_matrix <- cor_matrix
  }

  # Maximum correlation (should be with assigned cluster)
  result$max_cor <- round(apply(metric_matrix, 1, max), round_digits)

  # Second maximum correlation
  result$second_max_cor <- round(apply(metric_matrix, 1, function(x) {
    sorted <- sort(x, decreasing = TRUE)
    if (length(sorted) >= 2) sorted[2] else 0
  }), round_digits)

  # Separation: difference between max and second max
  # Higher values indicate clearer cluster assignment
  result$separation <- round(result$max_cor - result$second_max_cor, round_digits)

  # Reorder columns for clarity
  cluster_cols <- paste0("Cluster", 1:K)
  result <- result[, c("variable", "assigned_cluster", cluster_cols,
                       "max_cor", "second_max_cor", "separation")]

  # Sort by assigned cluster, then by separation (descending)
  result <- result[order(result$assigned_cluster, -result$separation), ]

  # Reset row names
  rownames(result) <- NULL

  return(result)
}
