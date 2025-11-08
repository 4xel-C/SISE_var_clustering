#' Visualization Functions for K-means Variable Clustering
#'
#' @description
#' Collection of plotting functions to visualize K-means variable clustering results.
#'
#' @details
#' This file provides:
#' - Elbow plot (inertia vs K)
#' - Silhouette plot
#' - Heatmap of correlation matrix (reorganized by clusters)
#' - Contribution barplots
#' - Diagnostic plots combining multiple metrics
#'
#' @note
#' These functions are internal utilities for the package.
#' They rely on base R graphics for minimal dependencies.
#'
#' @keywords internal
#' @noRd


#' Plot Elbow Curve
#'
#' @description
#' Visualizes the within-cluster inertia as a function of K.
#' Helps identify the "elbow point" where adding clusters provides diminishing returns.
#'
#' @param elbow_data data.frame from kmeans_elbow() with columns k, inertia
#' @param main Character. Plot title. Default: "Elbow Method"
#' @param xlab Character. X-axis label. Default: "Number of clusters (K)"
#' @param ylab Character. Y-axis label. Default: "Within-cluster inertia"
#' @param col Character. Line color. Default: "steelblue"
#' @param lwd Numeric. Line width. Default: 2
#' @param ... Additional arguments passed to plot()
#'
#' @return NULL (creates plot)
#'
#' @examples
#' \dontrun{
#' elbow_data <- kmeans_elbow(data, k_range = 2:10)
#' plot_kmeans_elbow(elbow_data)
#' }
#'
#' @noRd
plot_kmeans_elbow <- function(elbow_data, 
                              main = "Elbow Method for Optimal K",
                              xlab = "Number of clusters (K)",
                              ylab = "Within-cluster inertia",
                              col = "steelblue",
                              lwd = 2,
                              ...) {
  
  plot(elbow_data$k, elbow_data$inertia, 
       type = "b", 
       pch = 19,
       col = col,
       lwd = lwd,
       main = main,
       xlab = xlab,
       ylab = ylab,
       las = 1,
       ...)
  
  # Add grid for readability
  grid()
  
  # Add points again on top of grid
  points(elbow_data$k, elbow_data$inertia, pch = 19, col = col, cex = 1.2)
}


#' Plot Silhouette Scores with Variable Names
#'
#' @description
#' Visualizes silhouette coefficients for each variable, grouped by cluster.
#' Variables are sorted by silhouette score within each cluster.
#' Variable names are displayed on the y-axis for easy identification.
#'
#' @param silhouette_data data.frame from kmeans_silhouette() with columns:
#'        variable, cluster, silhouette
#' @param main Character. Plot title. Default: "Silhouette Plot"
#' @param col_palette Character vector. Colors for clusters.
#'        Default: rainbow palette
#' @param cex_names Numeric. Text size for variable names. Default: 0.7
#' @param ... Additional arguments passed to barplot()
#'
#' @return NULL (creates plot)
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' sil <- kmeans_silhouette(data, km$clusters, km$centroids)
#' plot_kmeans_silhouette(sil)
#' }
#'
#' @noRd
plot_kmeans_silhouette <- function(silhouette_data,
                                   main = "Silhouette Plot",
                                   col_palette = NULL,
                                   cex_names = 0.7,
                                   ...) {

  # Sort by cluster, then by silhouette value (descending)
  silhouette_data <- silhouette_data[order(silhouette_data$cluster,
                                           -silhouette_data$silhouette), ]

  # Color palette
  K <- max(silhouette_data$cluster)
  if (is.null(col_palette)) {
    col_palette <- rainbow(K, s = 0.6, v = 0.8)
  }

  # Assign colors
  colors <- col_palette[silhouette_data$cluster]

  # Adjust margins for variable names
  par(mar = c(5, 8, 4, 2) + 0.1)

  # Create horizontal barplot with variable names
  bp <- barplot(silhouette_data$silhouette,
          names.arg = silhouette_data$variable,
          col = colors,
          border = NA,
          main = main,
          xlab = "Silhouette coefficient",
          ylab = "",
          xlim = c(-1, 1),
          las = 1,
          horiz = TRUE,
          cex.names = cex_names,
          ...)

  # Add vertical line at 0
  abline(v = 0, lty = 2, col = "gray30", lwd = 1.5)

  # Add average silhouette line
  avg_sil <- mean(silhouette_data$silhouette)
  abline(v = avg_sil, lty = 2, col = "red", lwd = 2)

  # Add legend
  legend("bottomright",
         legend = c(paste0("Cluster ", 1:K),
                    "",
                    paste0("Moyenne: ", round(avg_sil, 3))),
         fill = c(col_palette, NA, NA),
         border = c(rep("black", K), NA, NA),
         lty = c(rep(NA, K), NA, 2),
         col = c(rep(NA, K), NA, "red"),
         lwd = c(rep(NA, K), NA, 2),
         bty = "n",
         cex = 0.8)
}


#' Plot Average Silhouette vs K
#'
#' @description
#' Visualizes average silhouette coefficient for different values of K.
#' Optimal K maximizes the average silhouette.
#'
#' @param silhouette_range_data data.frame from kmeans_silhouette_range()
#' @param main Character. Plot title.
#' @param col Character. Line color. Default: "darkgreen"
#' @param lwd Numeric. Line width. Default: 2
#' @param ... Additional arguments passed to plot()
#'
#' @return NULL (creates plot)
#'
#' @examples
#' \dontrun{
#' sil_data <- kmeans_silhouette_range(data, k_range = 2:8)
#' plot_kmeans_silhouette_range(sil_data)
#' }
#'
#' @noRd
plot_kmeans_silhouette_range <- function(silhouette_range_data,
                                         main = "Average Silhouette vs K",
                                         col = "darkgreen",
                                         lwd = 2,
                                         ...) {
  
  plot(silhouette_range_data$k, silhouette_range_data$avg_silhouette,
       type = "b",
       pch = 19,
       col = col,
       lwd = lwd,
       main = main,
       xlab = "Number of clusters (K)",
       ylab = "Average silhouette coefficient",
       las = 1,
       ylim = c(0, 1),
       ...)
  
  grid()
  points(silhouette_range_data$k, silhouette_range_data$avg_silhouette, 
         pch = 19, col = col, cex = 1.2)
  
  # Highlight optimal K
  optimal_k <- silhouette_range_data$k[which.max(silhouette_range_data$avg_silhouette)]
  optimal_sil <- max(silhouette_range_data$avg_silhouette)
  
  points(optimal_k, optimal_sil, pch = 19, col = "red", cex = 2)
  text(optimal_k, optimal_sil, 
       labels = paste0("K=", optimal_k), 
       pos = 3, col = "red", font = 2)
}


#' Plot Top Contributing Variables per Cluster
#'
#' @description
#' Visualizes the top contributing variables (highest R²) for each cluster.
#'
#' @param contribution_data data.frame from kmeans_contributions()
#' @param top_n Integer. Number of top variables to show per cluster. Default: 5
#' @param main Character. Plot title.
#' @param col_palette Character vector. Colors for clusters.
#' @param ... Additional arguments passed to barplot()
#'
#' @return NULL (creates plot)
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' contrib <- kmeans_contributions(data, km$clusters, km$centroids)
#' plot_kmeans_contributions(contrib, top_n = 5)
#' }
#'
#' @noRd
plot_kmeans_contributions <- function(contribution_data, 
                                      top_n = 5,
                                      main = "Top Contributing Variables per Cluster",
                                      col_palette = NULL,
                                      ...) {
  
  K <- max(contribution_data$cluster)
  
  # Color palette
  if (is.null(col_palette)) {
    col_palette <- rainbow(K, s = 0.6, v = 0.8)
  }
  
  # Setup multi-panel plot
  par(mfrow = c(1, K), mar = c(5, 8, 4, 2))
  
  for (k in 1:K) {
    
    # Filter and sort by contribution
    cluster_data <- contribution_data[contribution_data$cluster == k, ]
    cluster_data <- cluster_data[order(-cluster_data$contribution), ]
    
    # Select top N
    top_vars <- head(cluster_data, top_n)
    
    # Reverse order for barplot (top variable at top)
    top_vars <- top_vars[nrow(top_vars):1, ]
    
    # Create horizontal barplot
    barplot(top_vars$contribution,
            names.arg = top_vars$variable,
            horiz = TRUE,
            las = 1,
            col = col_palette[k],
            main = paste0("Cluster ", k),
            xlab = "Contribution (R²)",
            xlim = c(0, 1),
            ...)
  }
  
  # Reset graphical parameters
  par(mfrow = c(1, 1))
}


#' Create Diagnostic Plot for K-means Clustering
#'
#' @description
#' Creates a comprehensive 2x2 panel plot showing:
#' - Elbow curve (inertia vs K)
#' - Silhouette coefficient vs K
#' - Cluster sizes (if fitted_k provided)
#' - Intra-cluster correlation (if fitted_k provided)
#'
#' @param data data.frame or matrix of quantitative variables
#' @param k_range Integer vector of K values to test. Default: 2:10
#' @param fitted_k Integer. If provided, shows detailed stats for this K. Default: NULL
#' @param n_init Integer. Number of initializations. Default: 10
#' @param random_state Integer or NULL. Random seed.
#' @param ... Additional arguments
#'
#' @return NULL (creates plot)
#'
#' @examples
#' \dontrun{
#' plot_kmeans_diagnostics(data, k_range = 2:8, fitted_k = 5)
#' }
#'
#' @noRd
plot_kmeans_diagnostics <- function(data, k_range = 2:10, fitted_k = NULL,
                                    n_init = 10, random_state = NULL,
                                    correlation_type = "squared", ...) {

  # Compute metrics for range of K with specified correlation_type
  elbow_data <- kmeans_elbow(data, k_range, n_init, random_state,
                            correlation_type = correlation_type)
  sil_data <- kmeans_silhouette_range(data, k_range, n_init, random_state,
                                      correlation_type = correlation_type)

  # Setup 2x2 panel
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # Plot 1: Elbow
  plot(elbow_data$k, elbow_data$inertia, type = "b", pch = 19,
       col = "steelblue", lwd = 2,
       main = "Elbow Method",
       xlab = "K", ylab = "Inertia", las = 1)
  grid()
  if (!is.null(fitted_k)) {
    points(fitted_k, elbow_data$inertia[elbow_data$k == fitted_k],
           pch = 19, col = "red", cex = 2)
  }

  # Plot 2: Silhouette
  plot(sil_data$k, sil_data$avg_silhouette, type = "b", pch = 19,
       col = "darkgreen", lwd = 2,
       main = "Average Silhouette",
       xlab = "K", ylab = "Silhouette", las = 1, ylim = c(0, 1))
  grid()
  optimal_sil <- sil_data$k[which.max(sil_data$avg_silhouette)]
  points(optimal_sil, max(sil_data$avg_silhouette), pch = 19, col = "red", cex = 2)

  # Plot 3: Cluster sizes (if fitted_k provided)
  if (!is.null(fitted_k)) {
    km <- KmeansVariables$new(n_clusters = fitted_k, random_state = random_state,
                              correlation_type = correlation_type)
    km$fit(data)

    cluster_sizes <- as.numeric(table(km$clusters))
    barplot(cluster_sizes,
            names.arg = paste0("C", 1:fitted_k),
            col = rainbow(fitted_k, s = 0.6, v = 0.8),
            main = paste0("Cluster Sizes (K=", fitted_k, ")"),
            xlab = "Cluster", ylab = "N variables",
            las = 1)

    # Plot 4: Intra-cluster correlations
    intra_cor <- kmeans_intra_correlation(data, km$clusters)
    barplot(intra_cor,
            names.arg = paste0("C", 1:fitted_k),
            col = rainbow(fitted_k, s = 0.6, v = 0.8),
            main = "Intra-Cluster Correlation",
            xlab = "Cluster", ylab = "Mean |cor|",
            las = 1,
            ylim = c(0, 1))
    abline(h = 0.5, lty = 2, col = "red")

  } else {
    # Empty plots with messages
    plot.new()
    text(0.5, 0.5, "Specify 'fitted_k' to see\ncluster sizes", cex = 1.2)

    plot.new()
    text(0.5, 0.5, "Specify 'fitted_k' to see\nintra-cluster correlations", cex = 1.2)
  }

  # Reset graphical parameters
  par(mfrow = c(1, 1))
}


#' Plot Cluster Summary Panel
#'
#' @description
#' Creates a comprehensive visualization for a fitted K-means model showing:
#' - Silhouette plot with variable names
#' - Variable projection (PCA-based visualization)
#' - Top contributing variables per cluster
#'
#' @param km_model Fitted KmeansVariables object
#' @param data data.frame or matrix used for fitting
#' @param top_n Integer. Number of top contributors to show. Default: 5
#' @param ... Additional arguments
#'
#' @return NULL (creates plot)
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' plot_kmeans_summary(km, data)
#' }
#'
#' @noRd
plot_kmeans_summary <- function(km_model, data, top_n = 5, ...) {

  if (!km_model$fitted) {
    stop("Model must be fitted before plotting")
  }

  # Compute metrics
  sil <- kmeans_silhouette(data, km_model$clusters, km_model$centroids)
  contrib <- kmeans_contributions(data, km_model$clusters, km_model$centroids)

  # Setup layout: 2 rows, first row has 2 plots, second row spans full width
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))

  # Plot 1: Silhouette
  plot_kmeans_silhouette(sil, main = "Silhouette par Variable")

  # Plot 2: Variable Projection (PCA)
  plot_kmeans_projection(data, km_model$clusters,
                         main = "Projection des Variables (PCA)",
                         show_labels = TRUE)

  # Plot 3: Top Contributions
  plot_kmeans_contributions(contrib, top_n = top_n,
                           main = paste0("Top ", top_n, " Variables par Cluster"))

  # Reset layout
  layout(1)
}

#' Plot Variable Projection (PCA-style) 
#'
#' @noRd
plot_kmeans_projection <- function(data, clusters,
                                   main = "Variable Projection (PCA)",
                                   col_palette = NULL,
                                   show_labels = TRUE,
                                   label_threshold = 0.3,
                                   ...) {
  
  # PCA on correlation matrix
  cor_mat <- cor(data)
  pca <- prcomp(cor_mat, center = TRUE, scale. = TRUE)
  
  # Loadings scaled by sqrt(eigenvalues) = correlations with PCs
  # This gives values between -1 and +1
  loadings <- pca$rotation[, 1:2] %*% diag(sqrt(pca$sdev[1:2]^2))
  rownames(loadings) <- colnames(data)
  
  # Assign coords for easier reference
  coords <- loadings
  
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  # Quality of representation (cos²)
  # Distance to origin squared (on unit circle, 1 = perfect)
  quality <- rowSums(coords^2)
  
  # Color palette with better contrast
  K <- max(clusters)
  if (is.null(col_palette)) {
    col_palette <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6", "#1ABC9C")[1:K]
  }
  
  colors <- col_palette[clusters]
  
  # Point size based on quality (bounded between 0.5 and 2.5)
  cex_values <- pmax(0.5, pmin(2.5, 0.8 + 1.5 * quality))
  
  # PLOT WITH ASP = 1 (square aspect ratio)
  par(mar = c(5, 5, 4, 2), bg = "white")
  
  # Fixed limits for correlation circle: -1 to +1
  max_range <- 1.1
  
  plot(coords[, 1], coords[, 2],
       col = colors,
       pch = 19,
       cex = cex_values,
       xlab = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
       ylab = sprintf("PC2 (%.1f%% variance)", var_exp[2]),
       main = main,
       las = 1,
       xlim = c(-max_range, max_range),
       ylim = c(-max_range, max_range),
       asp = 1,  # Force aspect ratio = 1:1
       ...)
  
  # Grid
  abline(h = 0, v = 0, col = "gray40", lty = 2, lwd = 1.5)
  grid(col = "gray90")
  
  # Redraw points on top of grid
  points(coords[, 1], coords[, 2],
         col = colors,
         pch = 19,
         cex = cex_values)
  
  # Labels ONLY for very well-represented variables
  if (show_labels) {
    very_well_represented <- quality > label_threshold
    
    if (sum(very_well_represented) > 0) {
      text(coords[very_well_represented, 1], 
           coords[very_well_represented, 2],
           labels = rownames(coords)[very_well_represented],
           pos = 3,
           cex = 0.75,
           col = colors[very_well_represented],
           font = 2,
           offset = 0.5)
    }
  }
  
  # Correlation circle (radius = 1)
  theta <- seq(0, 2 * pi, length.out = 100)
  lines(cos(theta), sin(theta), 
        col = "gray30", lty = 3, lwd = 2)
  
  # Legend
  legend("topright",
         legend = c(paste0("Cluster ", 1:K), 
                    "",
                    sprintf("Quality threshold: %.2f", label_threshold),
                    sprintf("Labeled: %d/%d variables", 
                            sum(quality > label_threshold), nrow(coords))),
         fill = c(col_palette, rep(NA, 3)),
         border = c(rep("black", K), rep(NA, 3)),
         bty = "n",
         cex = 0.85)
  
  invisible(list(coords = coords, quality = quality, var_exp = var_exp))
}