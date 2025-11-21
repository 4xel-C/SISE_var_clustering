# ===========================================================================
# Visualization Functions for K-means Variable Clustering
# ===========================================================================
#
# Essential plotting functions to visualize K-means variable clustering:
# - plot_kmeans_elbow: Elbow plot (homogeneity vs K) for choosing optimal K
# - plot_kmeans_projection: PCA projection plot for visualizing variable groups
# - plot_kmeans_contributions: Top contributing variables per cluster
#
# All functions are internal utilities (@noRd) and rely on base R graphics.
# ===========================================================================


#' Plot Elbow Curve
#'
#' @description
#' Visualizes the within-cluster homogeneity (cohesion) as a function of K.
#' Helps identify the "elbow point" where adding clusters provides diminishing
#' returns.
#'
#' @param elbow_data data.frame from kmeans_elbow() with columns k, inertia
#' @param main Character. Plot title. Default: "Elbow Method for Optimal K"
#' @param xlab Character. X-axis label. Default: "Number of clusters (K)"
#' @param ylab Character. Y-axis label. Default: "Homogeneity criterion"
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
                              ylab = "Homogeneity criterion",
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
       xaxt = "n",   # Désactive l’axe x
       ...
  )

  axis(1, at = elbow_data$k, labels = elbow_data$k, las = 1)

  # Add grid for readability
  grid()

  # Add points again on top of grid
  points(elbow_data$k, elbow_data$inertia, pch = 19, col = col, cex = 1.2)
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
    top_vars <- top_vars[rev(seq_len(nrow(top_vars))), ]

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


#' Plot Variable Projection (PCA-style)
#'
#' @description
#' Projects variables onto the first two principal components of the correlation
#' matrix. Variables are colored by cluster assignment. This visualization helps
#' understand the geometric structure of clusters and identify well-represented
#' variables.
#'
#' @param data data.frame or matrix of quantitative variables
#' @param clusters Integer vector. Cluster assignment for each variable
#' @param main Character. Plot title. Default: "Variable Projection (PCA)"
#' @param col_palette Character vector. Colors for clusters. Default: NULL
#' @param show_labels Logical. Show variable names? Default: TRUE
#' @param label_threshold Numeric. Only label variables with quality > threshold.
#'        Set to 0 to display all labels. Default: 0
#' @param supplementary_variables Optional data.frame/matrix of supplementary variables
#'        to project onto the PCA space. Default: NULL
#' @param supplementary_clusters Optional integer vector of cluster assignments
#'        for supplementary variables (from predict()). Default: NULL
#' @param ... Additional arguments passed to plot()
#'
#' @return Invisibly returns a list with coords, quality, and var_exp
#'
#' @examples
#' \dontrun{
#' km <- KmeansVariables$new(n_clusters = 3)
#' km$fit(data)
#' plot_kmeans_projection(data, km$clusters, show_labels = TRUE)
#' }
#'
#' @noRd
plot_kmeans_projection <- function(data, clusters,
                                   main = "Variable Projection (PCA)",
                                   col_palette = NULL,
                                   show_labels = TRUE,
                                   label_threshold = 0,
                                   supplementary_variables = NULL,
                                   supplementary_clusters = NULL,
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
    col_palette <- c("#E74C3C", "#3498DB", "#2ECC71",
                     "#F39C12", "#9B59B6", "#1ABC9C")[1:K]
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

  # ===========================================================================
  # PROJECT SUPPLEMENTARY VARIABLES (if provided)
  # ===========================================================================

  if (!is.null(supplementary_variables)) {

    X_sup <- as.matrix(supplementary_variables)

    # Check same number of observations as training data
    if (nrow(X_sup) != nrow(data)) {
      warning("Supplementary variables must have same number of observations. Skipping projection.")
    } else {

      # Calculate correlations: supplementary vars × training vars
      cor_sup_train <- cor(X_sup, as.matrix(data))

      # Project onto existing PCA axes
      coords_sup <- cor_sup_train %*% pca$rotation[, 1:2] %*% diag(sqrt(pca$sdev[1:2]^2))
      rownames(coords_sup) <- colnames(X_sup)

      # Quality of representation
      quality_sup <- rowSums(coords_sup^2)

      # Colors: use cluster colors if provided, else gray
      if (!is.null(supplementary_clusters)) {
        colors_sup <- col_palette[supplementary_clusters]
      } else {
        colors_sup <- rep("darkgray", nrow(coords_sup))
      }

      # Point size based on quality
      cex_sup <- pmax(0.5, pmin(2.5, 0.8 + 1.5 * quality_sup))

      # Plot with TRIANGLES (to distinguish from training vars)
      points(coords_sup[, 1], coords_sup[, 2],
             col = colors_sup,
             pch = 17,  # Triangle symbol
             cex = cex_sup)

      # Labels for well-represented supplementary variables
      if (show_labels) {
        well_represented_sup <- quality_sup > label_threshold

        if (sum(well_represented_sup) > 0) {
          text(coords_sup[well_represented_sup, 1],
               coords_sup[well_represented_sup, 2],
               labels = rownames(coords_sup)[well_represented_sup],
               pos = 1,  # Below point
               cex = 0.7,
               col = colors_sup[well_represented_sup],
               font = 3)  # Italic font
        }
      }
    }
  }

  # Correlation circle (radius = 1)
  theta <- seq(0, 2 * pi, length.out = 100)
  lines(cos(theta), sin(theta),
        col = "gray30", lty = 3, lwd = 2)

  # ===========================================================================
  # LEGEND (CORRECTED - two separate legends)
  # ===========================================================================

  # Legend 1: Clusters with colors
  legend("topright",
         legend = paste0("Cluster ", 1:K),
         fill = col_palette,
         border = "black",
         bty = "n",
         cex = 0.85)

  # Legend 2: Info text (positioned below first legend)
  n_labeled <- sum(quality > label_threshold)
  legend("topright",
         legend = c("",
                    sprintf("Quality threshold: %.2f", label_threshold),
                    sprintf("Labeled: %d/%d variables", n_labeled, nrow(coords))),
         bty = "n",
         cex = 0.75,
         text.col = "gray40",
         inset = c(0, 0.05 * (K + 1)))  # Adjust position based on K

  invisible(list(coords = coords, quality = quality, var_exp = var_exp))
}
