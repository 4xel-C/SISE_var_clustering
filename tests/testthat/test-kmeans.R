#' Unit Tests for K-means Variable Clustering
#'
#' @description
#' Test suite for the KmeansVariables class and associated functions.
#' Uses testthat framework for comprehensive testing.
#'
#' @details
#' Tests cover:
#' - Object initialization
#' - Parameter validation
#' - Fitting on synthetic data
#' - Prediction on new variables
#' - Edge cases and error handling
#' - Metrics calculation
#' - Reproducibility with random_state
#'
#' @note
#' These tests should be placed in tests/testthat/ directory
#' and named test-kmeans.R to be automatically discovered by testthat.
#'
#' @keywords internal


# =============================================================================
# TEST SETUP
# =============================================================================

# Helper function to create synthetic clustered data
create_synthetic_data <- function(n_vars_per_cluster = 5, 
                                  n_obs = 500, 
                                  n_clusters = 3,
                                  cor_strengths = c(0.7, 0.65, 0.6)) {
  
  if (length(cor_strengths) != n_clusters) {
    stop("cor_strengths must have length equal to n_clusters")
  }
  
  set.seed(123)
  data_list <- list()
  
  for (k in 1:n_clusters) {
    # Create covariance matrix
    cov_matrix <- matrix(cor_strengths[k], n_vars_per_cluster, n_vars_per_cluster)
    diag(cov_matrix) <- 1
    
    # Generate data
    cluster_data <- MASS::mvrnorm(n_obs, rep(0, n_vars_per_cluster), cov_matrix)
    colnames(cluster_data) <- paste0("Var", k, "_", 1:n_vars_per_cluster)
    data_list[[k]] <- cluster_data
  }
  
  # Combine into single data frame
  data <- as.data.frame(do.call(cbind, data_list))
  
  # True labels for validation
  true_labels <- rep(1:n_clusters, each = n_vars_per_cluster)
  
  return(list(data = data, true_labels = true_labels))
}


# =============================================================================
# INITIALIZATION TESTS
# =============================================================================

test_that("KmeansVariables initializes with default parameters", {
  
  km <- KmeansVariables$new()
  
  expect_s3_class(km, "KmeansVariables")
  expect_s3_class(km, "ClusteringBase")
  expect_equal(km$n_clusters, 3)
  expect_false(km$fitted)
})


test_that("KmeansVariables initializes with custom parameters", {
  
  km <- KmeansVariables$new(
    n_clusters = 5,
    max_iter = 200,
    tol = 1e-5,
    n_init = 20,
    random_state = 42
  )
  
  expect_equal(km$n_clusters, 5)
  expect_false(km$fitted)
})


test_that("KmeansVariables rejects invalid parameters", {
  
  # Invalid n_clusters
  expect_error(KmeansVariables$new(n_clusters = 0))
  expect_error(KmeansVariables$new(n_clusters = -1))
  expect_error(KmeansVariables$new(n_clusters = "three"))
  
  # Invalid max_iter
  expect_error(KmeansVariables$new(max_iter = 0))
  expect_error(KmeansVariables$new(max_iter = -10))
  
  # Invalid tol
  expect_error(KmeansVariables$new(tol = 0))
  expect_error(KmeansVariables$new(tol = -1e-4))
  
  # Invalid n_init
  expect_error(KmeansVariables$new(n_init = 0))
  expect_error(KmeansVariables$new(n_init = -5))
})


# =============================================================================
# FIT METHOD TESTS
# =============================================================================

test_that("fit() works on synthetic data with known clusters", {
  
  # Create synthetic data with 3 clusters
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  true_labels <- synthetic$true_labels
  
  # Fit model
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  # Check fitted status
  expect_true(km$fitted)
  
  # Check outputs exist
  expect_length(km$clusters, ncol(data))
  expect_equal(nrow(km$centroids), nrow(data))
  expect_equal(ncol(km$centroids), 3)
  expect_type(km$inertia, "double")
  expect_gte(km$n_iter, 1)
  
  # Check cluster sizes
  expect_equal(length(unique(km$clusters)), 3)
  
  # Check inertia is reasonable (should be relatively low for well-separated clusters)
  expect_lt(km$inertia, ncol(data) * 0.5)  # Less than 50% of max possible inertia
})


test_that("fit() produces consistent results with random_state", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  # Fit twice with same seed
  km1 <- KmeansVariables$new(n_clusters = 3, random_state = 42)
  km1$fit(data)
  
  km2 <- KmeansVariables$new(n_clusters = 3, random_state = 42)
  km2$fit(data)
  
  # Results should be identical
  expect_equal(km1$clusters, km2$clusters)
  expect_equal(km1$inertia, km2$inertia)
  expect_equal(km1$n_iter, km2$n_iter)
})


test_that("fit() rejects invalid data", {
  
  km <- KmeansVariables$new(n_clusters = 3)
  
  # Non-dataframe/matrix
  expect_error(km$fit("not a dataframe"))
  expect_error(km$fit(list(a = 1, b = 2)))
  
  # Too few variables
  expect_error(km$fit(data.frame(x = 1:10)))
  
  # Data with missing values
  data_na <- data.frame(x = c(1, 2, NA, 4), y = c(5, 6, 7, 8))
  expect_error(km$fit(data_na))
  
  # Non-numeric data
  data_char <- data.frame(x = letters[1:10], y = letters[11:20])
  expect_error(km$fit(data_char))
})


test_that("fit() rejects n_clusters > n_variables", {
  
  data <- data.frame(x = 1:10, y = 11:20, z = 21:30)
  
  km <- KmeansVariables$new(n_clusters = 5)  # More clusters than variables
  expect_error(km$fit(data))
})


test_that("fit() handles single variable per cluster", {
  
  # Create data where K = p (one variable per cluster)
  data <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100),
    x3 = rnorm(100)
  )
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  expect_no_error(km$fit(data))
  expect_true(km$fitted)
})


test_that("fit() converges within max_iter", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, max_iter = 50, random_state = 123)
  km$fit(data)
  
  # Should converge before max_iter for well-separated clusters
  expect_lte(km$n_iter, 50)
})


# =============================================================================
# PREDICT METHOD TESTS
# =============================================================================

test_that("predict() assigns new variables to clusters", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  # Fit model
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  # Create new variables (similar to cluster 1)
  new_data <- data.frame(
    NewVar1 = data[, 1] + rnorm(nrow(data), 0, 0.1),
    NewVar2 = data[, 2] + rnorm(nrow(data), 0, 0.1)
  )
  
  # Predict
  predictions <- km$predict(new_data)
  
  # Check output structure
  expect_type(predictions, "list")
  expect_named(predictions, c("clusters", "distances", "correlations"))
  expect_length(predictions$clusters, 2)
  expect_equal(nrow(predictions$distances), 2)
  expect_equal(ncol(predictions$distances), 3)
})


test_that("predict() requires fitted model", {
  
  km <- KmeansVariables$new(n_clusters = 3)
  new_data <- data.frame(x = 1:10, y = 11:20)
  
  expect_error(km$predict(new_data), "Model must be fitted")
})


test_that("predict() rejects mismatched number of observations", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  # Wrong number of observations
  new_data_wrong <- data.frame(x = 1:100, y = 1:100)  # 100 obs instead of 500
  
  expect_error(km$predict(new_data_wrong), "must have")
})


# =============================================================================
# PRINT AND SUMMARY TESTS
# =============================================================================

test_that("print() works for unfitted model", {
  
  km <- KmeansVariables$new(n_clusters = 3)
  
  expect_output(km$print(), "not fitted")
})


test_that("print() displays model information", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  expect_output(km$print(), "K-means Variable Clustering")
  expect_output(km$print(), "Number of clusters")
  expect_output(km$print(), "Iterations")
  expect_output(km$print(), "Total inertia")
})


test_that("summary() provides detailed information", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  expect_output(km$summary(), "Cluster Statistics")
  expect_output(km$summary(), "Within-cluster inertia")
  expect_output(km$summary(), "Mean correlation")
  expect_output(km$summary(), "Top variables")
})


# =============================================================================
# ACTIVE BINDINGS TESTS
# =============================================================================

test_that("Active bindings return NULL for unfitted model", {
  
  km <- KmeansVariables$new(n_clusters = 3)
  
  expect_warning(km$clusters)
  expect_warning(km$centroids)
  expect_warning(km$inertia)
  expect_warning(km$n_iter)
  expect_warning(km$cluster_sizes)
  expect_warning(km$cluster_inertias)
})


test_that("Active bindings return correct values after fitting", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  # clusters
  expect_length(km$clusters, 15)
  expect_true(all(km$clusters %in% 1:3))
  
  # centroids
  expect_equal(dim(km$centroids), c(500, 3))
  
  # inertia
  expect_type(km$inertia, "double")
  expect_gte(km$inertia, 0)
  
  # n_iter
  expect_type(km$n_iter, "integer")
  expect_gte(km$n_iter, 1)
  
  # cluster_sizes
  expect_length(km$cluster_sizes, 3)
  expect_equal(sum(km$cluster_sizes), 15)
  
  # cluster_inertias
  expect_length(km$cluster_inertias, 3)
  expect_true(all(km$cluster_inertias >= 0))
})


# =============================================================================
# METRICS TESTS
# =============================================================================

test_that("kmeans_elbow() computes inertia for range of K", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  elbow_data <- kmeans_elbow(data, k_range = 2:5, random_state = 123)
  
  expect_s3_class(elbow_data, "data.frame")
  expect_named(elbow_data, c("k", "inertia", "inertia_pct"))
  expect_equal(nrow(elbow_data), 4)
  
  # Inertia should decrease as K increases
  expect_true(all(diff(elbow_data$inertia) < 0))
})


test_that("kmeans_silhouette() computes silhouette coefficients", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  sil <- kmeans_silhouette(data, km$clusters, km$centroids)
  
  expect_s3_class(sil, "data.frame")
  expect_named(sil, c("variable", "cluster", "silhouette"))
  expect_equal(nrow(sil), ncol(data))
  
  # Silhouette should be between -1 and 1
  expect_true(all(sil$silhouette >= -1 & sil$silhouette <= 1))
  
  # For well-separated clusters, average silhouette should be positive
  expect_gt(mean(sil$silhouette), 0)
})


test_that("kmeans_calinski_harabasz() computes CH index", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  ch <- kmeans_calinski_harabasz(data, km$clusters, km$centroids)
  
  expect_type(ch, "double")
  expect_gt(ch, 0)  # CH should be positive
})


test_that("kmeans_contributions() computes variable contributions", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  contrib <- kmeans_contributions(data, km$clusters, km$centroids)
  
  expect_s3_class(contrib, "data.frame")
  expect_named(contrib, c("variable", "cluster", "contribution", "correlation"))
  expect_equal(nrow(contrib), ncol(data))
  
  # Contributions should be between 0 and 1
  expect_true(all(contrib$contribution >= 0 & contrib$contribution <= 1))
})


test_that("kmeans_intra_correlation() computes within-cluster correlations", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  intra_cor <- kmeans_intra_correlation(data, km$clusters)
  
  expect_length(intra_cor, 3)
  expect_named(intra_cor, c("Cluster1", "Cluster2", "Cluster3"))
  
  # Correlations should be between 0 and 1 (absolute values)
  expect_true(all(intra_cor >= 0 & intra_cor <= 1, na.rm = TRUE))
})


test_that("kmeans_find_optimal_k() suggests reasonable K", {

  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data

  result <- kmeans_find_optimal_k(data, k_range = 2:5, random_state = 123)

  expect_type(result, "list")
  expect_named(result, c("optimal_k", "metrics", "method"))

  # Optimal K should be within the tested range
  expect_true(result$optimal_k %in% 2:5)
  expect_gte(result$optimal_k, 2)
  expect_lte(result$optimal_k, 5)

  # Metrics should be a data frame
  expect_s3_class(result$metrics, "data.frame")
  expect_true("k" %in% names(result$metrics))
  expect_true("inertia" %in% names(result$metrics))
  expect_true("avg_silhouette" %in% names(result$metrics))
})


# =============================================================================
# EDGE CASES AND ROBUSTNESS TESTS
# =============================================================================

test_that("Algorithm handles variables with perfect correlation", {
  
  # Create data where two variables are identical
  data <- data.frame(
    x1 = rnorm(100),
    x2 = rnorm(100),
    x3 = rnorm(100)
  )
  data$x4 <- data$x1  # Perfect correlation
  
  km <- KmeansVariables$new(n_clusters = 2, random_state = 123)
  expect_no_error(km$fit(data))
  
  # x1 and x4 should be in the same cluster
  expect_equal(km$clusters[1], km$clusters[4])
})


test_that("Algorithm handles variables with zero variance", {

  # Create data with one constant variable
  data <- data.frame(
    x1 = rnorm(100),
    x2 = rep(5, 100),  # Zero variance
    x3 = rnorm(100)
  )

  km <- KmeansVariables$new(n_clusters = 2, random_state = 123)

  # Should handle gracefully without errors
  expect_no_error(km$fit(data))
  expect_true(km$fitted)
})


test_that("Multiple initializations improve results", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  # Single initialization
  km1 <- KmeansVariables$new(n_clusters = 3, n_init = 1, random_state = 123)
  km1$fit(data)
  
  # Multiple initializations
  km10 <- KmeansVariables$new(n_clusters = 3, n_init = 10, random_state = 123)
  km10$fit(data)
  
  # Multiple inits should find better or equal solution
  expect_lte(km10$inertia, km1$inertia)
})


test_that("Algorithm handles large number of variables", {
  
  # Create dataset with many variables
  set.seed(123)
  data <- as.data.frame(matrix(rnorm(100 * 50), nrow = 100, ncol = 50))
  
  km <- KmeansVariables$new(n_clusters = 5, random_state = 123)
  expect_no_error(km$fit(data))
  expect_true(km$fitted)
})


# =============================================================================
# INTEGRATION TESTS
# =============================================================================

test_that("Full workflow: fit, predict, visualize", {
  
  synthetic <- create_synthetic_data(n_vars_per_cluster = 5, n_clusters = 3)
  data <- synthetic$data
  
  # Fit
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  # Print and summary
  expect_output(km$print())
  expect_output(km$summary())
  
  # Metrics
  sil <- kmeans_silhouette(data, km$clusters, km$centroids)
  expect_gt(mean(sil$silhouette), 0)
  
  ch <- kmeans_calinski_harabasz(data, km$clusters, km$centroids)
  expect_gt(ch, 0)
  
  contrib <- kmeans_contributions(data, km$clusters, km$centroids)
  expect_equal(nrow(contrib), ncol(data))
  
  # Predict
  new_data <- data[, 1:2] + matrix(rnorm(nrow(data) * 2, 0, 0.1), ncol = 2)
  predictions <- km$predict(new_data)
  expect_length(predictions$clusters, 2)
})


test_that("Clustering recovers true structure in synthetic data", {
  
  # Create highly separated clusters
  synthetic <- create_synthetic_data(
    n_vars_per_cluster = 5,
    n_clusters = 3,
    cor_strengths = c(0.9, 0.85, 0.8)
  )
  data <- synthetic$data
  true_labels <- synthetic$true_labels
  
  # Fit model
  km <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km$fit(data)
  
  # Calculate adjusted Rand index (requires external package, so simplified check)
  # For well-separated clusters, most variables should be correctly clustered
  # We check if variables from the same true cluster tend to be in the same predicted cluster
  
  cluster_purity <- numeric(3)
  for (k in 1:3) {
    true_cluster_vars <- which(true_labels == k)
    predicted_clusters <- km$clusters[true_cluster_vars]
    # Most common predicted cluster
    cluster_purity[k] <- max(table(predicted_clusters)) / length(predicted_clusters)
  }
  
  # Average purity should be high (>0.7) for well-separated clusters
  expect_gt(mean(cluster_purity), 0.7)
})


# =============================================================================
# AUTOMATIC K SELECTION TESTS
# =============================================================================

test_that("n_clusters = 'auto' works with default parameters", {

  # Create synthetic data with clear structure
  synthetic <- create_synthetic_data(
    n_vars_per_cluster = 5,
    n_clusters = 3,
    cor_strengths = c(0.8, 0.75, 0.7)
  )
  data <- synthetic$data

  # Fit with automatic K selection
  km_auto <- KmeansVariables$new(n_clusters = "auto", random_state = 123)

  # Should not throw error
  expect_no_error(km_auto$fit(data))

  # Should have selected a valid K
  expect_true(km_auto$fitted)
  expect_true(is.numeric(km_auto$n_clusters))
  expect_gte(km_auto$n_clusters, 2)
  expect_lte(km_auto$n_clusters, 10)

  # Should have valid results
  expect_length(km_auto$clusters, ncol(data))
  expect_true(!is.null(km_auto$centroids))
  expect_true(!is.null(km_auto$inertia))
})


test_that("n_clusters = 'auto' respects k_range parameter", {

  # Create small dataset
  data_small <- create_synthetic_data(n_vars_per_cluster = 3, n_clusters = 2)$data

  # Specify custom k_range
  km_auto <- KmeansVariables$new(
    n_clusters = "auto",
    k_range = 2:4,
    random_state = 123
  )
  km_auto$fit(data_small)

  # Selected K should be within specified range
  expect_gte(km_auto$n_clusters, 2)
  expect_lte(km_auto$n_clusters, 4)
})


test_that("n_clusters = 'auto' works with different selection methods", {

  data <- create_synthetic_data(n_vars_per_cluster = 4, n_clusters = 3)$data

  # Test with silhouette method
  km_sil <- KmeansVariables$new(
    n_clusters = "auto",
    k_range = 2:5,
    selection_method = "silhouette",
    random_state = 123
  )
  expect_no_error(km_sil$fit(data))
  expect_true(km_sil$fitted)

  # Test with calinski method
  km_ch <- KmeansVariables$new(
    n_clusters = "auto",
    k_range = 2:5,
    selection_method = "calinski",
    random_state = 123
  )
  expect_no_error(km_ch$fit(data))
  expect_true(km_ch$fitted)

  # Test with "all" method
  km_all <- KmeansVariables$new(
    n_clusters = "auto",
    k_range = 2:5,
    selection_method = "all",
    random_state = 123
  )
  expect_no_error(km_all$fit(data))
  expect_true(km_all$fitted)
})


test_that("n_clusters = 'auto' validates parameters correctly", {

  # Invalid k_range (less than 2 values)
  expect_error(
    KmeansVariables$new(n_clusters = "auto", k_range = 3),
    "k_range.*at least 2 values"
  )

  # Invalid k_range (contains values < 2)
  expect_error(
    KmeansVariables$new(n_clusters = "auto", k_range = c(1, 2, 3)),
    "k_range.*>= 2"
  )

  # Invalid selection_method
  expect_error(
    KmeansVariables$new(n_clusters = "auto", selection_method = "invalid"),
    "selection_method.*silhouette.*calinski.*all"
  )
})


test_that("Manual K selection still works as before", {

  # Ensure backward compatibility: manual K selection unchanged
  data <- create_synthetic_data(n_vars_per_cluster = 4, n_clusters = 3)$data

  km_manual <- KmeansVariables$new(n_clusters = 3, random_state = 123)
  km_manual$fit(data)

  # Should work exactly as before
  expect_equal(km_manual$n_clusters, 3)
  expect_true(km_manual$fitted)
  expect_length(km_manual$clusters, ncol(data))
})