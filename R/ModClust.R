#' Modality-Based Hierarchical Clustering with Dice Distance
#'
#' @description
#' \code{ModCluster} implements a modality-based clustering approach using the Dice
#' distance coefficient for categorical data. This method clusters modalities (categories
#' or levels of categorical variables) rather than individual observations, making it
#' particularly suitable for analyzing patterns in categorical data structures.
#'
#' The algorithm performs hierarchical clustering on a complete disjunctive table
#' (indicator matrix) using the Dice distance between modalities. It integrates
#' Multiple Correspondence Analysis (MCA) for visualization purposes and provides
#' prediction capabilities for new observations based on their modality profiles.
#'
#' @section Key Features:
#' \itemize{
#'   \item \strong{Dice Distance}: Uses the Dice dissimilarity coefficient specifically
#'         designed for binary and categorical data comparison
#'   \item \strong{Hierarchical Clustering}: Supports various linkage methods including
#'         Ward's method, average linkage, single linkage, complete linkage, and others
#'   \item \strong{MCA Visualization}: Integrates Multiple Correspondence Analysis to
#'         provide graphical representations of modality clusters in reduced dimensional space
#'   \item \strong{Prediction}: Predicts cluster membership for new observations based on
#'         their distance to existing modality clusters
#'   \item \strong{Projection}: Projects new observations into the MCA factorial space
#'         for visual interpretation
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(n_clusters, method, n_dimensions, hclust_method)}}{
#'     Initialize a new ModCluster object.
#'     \itemize{
#'       \item \code{n_clusters}: Integer, target number of clusters (default: NULL).
#'             If NULL, the number of clusters must be specified later using \code{cut_tree(k)}.
#'       \item \code{method}: Character, clustering method. Currently only "hclust" is supported
#'             (default: "hclust").
#'       \item \code{n_dimensions}: Integer, number of MCA dimensions to compute for visualization
#'             purposes (default: 5). Higher values provide more detailed representation but
#'             increase computational cost.
#'       \item \code{hclust_method}: Character, linkage method for hierarchical clustering.
#'             Must be one of: "ward.D2" (Ward's minimum variance), "ward.D", "single" (nearest neighbor),
#'             "complete" (farthest neighbor), "average" (UPGMA), "mcquitty" (WPGMA),
#'             "median" (WPGMC), "centroid" (UPGMC). Default is "average".
#'     }
#'   }
#'
#'   \item{\code{fit(data)}}{
#'     Fit the clustering model on the training data.
#'     \itemize{
#'       \item \code{data}: Data frame containing categorical variables. All variables
#'             should be factors or will be automatically converted to factors. The method
#'             performs the following steps: (1) creates a complete disjunctive table,
#'             (2) computes Dice distances between all modalities, (3) applies hierarchical
#'             clustering, and (4) computes MCA for visualization.
#'     }
#'     Returns: Self (invisibly) for method chaining.
#'   }
#'
#'   \item{\code{predict(new_data)}}{
#'     Predict cluster membership for new observations based on their modality profiles.
#'     \itemize{
#'       \item \code{new_data}: Data frame with the same categorical variables as the
#'             training data. Factor levels must match those in the training set.
#'     }
#'     Returns: Data frame with the following columns:
#'     \itemize{
#'       \item \code{observation}: Observation index (1 to n)
#'       \item \code{predicted_cluster}: Predicted cluster assignment (integer from 1 to k)
#'       \item \code{dist_cluster_1}, \code{dist_cluster_2}, ...: Dice distance to each
#'             cluster centroid, providing confidence information about the assignment
#'     }
#'   }
#'
#'   \item{\code{project_new_modalities(new_data, axes)}}{
#'     Project new observations into the MCA factorial space.
#'     \itemize{
#'       \item \code{new_data}: Data frame with the same categorical variables as training data.
#'       \item \code{axes}: Integer vector specifying which MCA dimensions to use for projection
#'             (default: all available dimensions).
#'     }
#'     Returns: Matrix of projected coordinates with rows representing observations and
#'     columns representing MCA dimensions (named "Dim.1", "Dim.2", etc.).
#'   }
#'
#'   \item{\code{cut_tree(k)}}{
#'     Cut the hierarchical tree into k clusters. This method must be called if n_clusters
#'     was not specified during initialization.
#'     \itemize{
#'       \item \code{k}: Integer, desired number of clusters.
#'     }
#'     Returns: Self (invisibly) for method chaining.
#'   }
#'
#'   \item{\code{plot_dendrogram(k, color_branches)}}{
#'     Plot the hierarchical clustering dendrogram showing the tree structure.
#'     \itemize{
#'       \item \code{k}: Integer, number of clusters to highlight with rectangles or colors
#'             (default: uses current k value if available).
#'       \item \code{color_branches}: Logical, whether to color dendrogram branches by cluster
#'             (default: TRUE). Requires the \pkg{dendextend} package for colored output.
#'     }
#'   }
#'
#'   \item{\code{plot_mca(axes, show_labels, label_cex)}}{
#'     Plot modalities in the MCA factorial space with colors representing cluster membership.
#'     \itemize{
#'       \item \code{axes}: Integer vector of length 2, specifying which MCA axes to plot
#'             (default: c(1, 2) for the first factorial plane).
#'       \item \code{show_labels}: Logical, whether to display modality labels on the plot
#'             (default: TRUE).
#'       \item \code{label_cex}: Numeric, character expansion factor for label text size
#'             (default: 0.7).
#'     }
#'   }
#'
#'   \item{\code{plot_mca_with_new(new_data, axes, show_labels, label_cex)}}{
#'     Plot both training modalities (as points) and new observations (as stars) in the
#'     MCA factorial space.
#'     \itemize{
#'       \item \code{new_data}: Data frame containing new observations to project and display.
#'       \item \code{axes}: Integer vector of length 2, specifying which MCA axes to plot
#'             (default: c(1, 2)).
#'       \item \code{show_labels}: Logical, whether to show labels for both training and
#'             new observations (default: TRUE).
#'       \item \code{label_cex}: Numeric, character expansion factor for label text
#'             (default: 0.7).
#'     }
#'     Returns: Invisible list containing:
#'     \itemize{
#'       \item \code{train_coords}: Matrix of training modality coordinates
#'       \item \code{new_coords}: Matrix of new observation coordinates
#'       \item \code{train_clusters}: Cluster assignments for training modalities
#'       \item \code{new_clusters}: Predicted cluster assignments for new observations
#'     }
#'   }
#'
#'   \item{\code{plot_silhouette()}}{
#'     Plot silhouette values for each modality to assess clustering quality. The silhouette
#'     coefficient measures how similar a modality is to its own cluster compared to other
#'     clusters. Values range from -1 to 1, where high values indicate good cluster assignment.
#'   }
#'
#'   \item{\code{plot_heights()}}{
#'     Plot aggregation heights from the hierarchical clustering to help determine the
#'     optimal number of clusters using the elbow method. The plot shows how the within-cluster
#'     dissimilarity changes as clusters are merged.
#'   }
#'
#'   \item{\code{summary()}}{
#'     Get a comprehensive summary of clustering results including cluster statistics and
#'     modality-level details.
#'     Returns: List with two data frames:
#'     \itemize{
#'       \item \code{clust_summary}: Cluster-level summary containing cluster ID and number
#'             of member modalities
#'       \item \code{clust_members}: Modality-level details including modality name, cluster
#'             assignment, frequency in the data, distance to own cluster centroid, distance
#'             to nearest other cluster, and the ratio of these distances
#'     }
#'   }
#'
#'   \item{\code{get_cluster_table()}}{
#'     Get a simple table of modalities with their cluster assignments and frequencies.
#'     Returns: Data frame with columns Modality (character), Cluster (integer), and
#'     Frequency (numeric).
#'   }
#'
#'   \item{\code{get_dice_matrix()}}{
#'     Get the complete Dice distance matrix between all modalities.
#'     Returns: Symmetric matrix of Dice distances with rows and columns labeled by
#'     modality names.
#'   }
#'
#'   \item{\code{print()}}{
#'     Print a concise summary of the ModCluster object showing the clustering method,
#'     number of modalities, number of clusters, and cluster sizes.
#'   }
#' }
#'
#' @section Active Bindings:
#' These read-only properties provide access to internal clustering results:
#' \describe{
#'   \item{\code{modality_labels}}{Named integer vector of cluster assignments for each modality.
#'         Names correspond to modality names in "variable.level" format.}
#'   \item{\code{modality_names}}{Character vector of all modality names in the format
#'         "variable.level" (e.g., "color.red", "size.large").}
#'   \item{\code{modality_frequencies}}{Named numeric vector showing the frequency (count)
#'         of each modality in the training data.}
#'   \item{\code{mca_result}}{Complete MCA result object from FactoMineR::MCA, containing
#'         eigenvalues, coordinates, contributions, and other MCA outputs.}
#'   \item{\code{cluster_centers}}{Matrix of cluster centers represented as binary profiles
#'         in the complete disjunctive table space.}
#' }
#'
#' @section Algorithm Details:
#' The Dice distance between two modalities j and j' is computed as:
#' \deqn{D(j, j') = \sqrt{\frac{1}{2} \sum_{i=1}^{n} (x_{ij}(1-x_{ij'}) + (1-x_{ij})x_{ij'})}}
#'
#' where \eqn{x_{ij}} is the indicator (0 or 1) that individual i has modality j, and n
#' is the total number of individuals.
#'
#' This distance measures the dissimilarity between two modalities based on how often they
#' co-occur or occur separately across individuals. It is particularly well-suited for
#' binary data and categorical variables represented as indicator matrices.
#'
#' The prediction for new observations is based on computing the average Dice distance
#' (or minimum/maximum depending on the linkage method) from the observation to all
#' modalities within each cluster, then assigning the observation to the cluster with
#' the smallest distance.
#'
#' @section Workflow:
#' A typical workflow consists of:
#' \enumerate{
#'   \item Initialize the model: \code{model <- ModCluster$new(n_clusters = 3)}
#'   \item Fit to training data: \code{model$fit(data)}
#'   \item Visualize results: \code{model$plot_dendrogram()}, \code{model$plot_mca()}
#'   \item Evaluate quality: \code{model$plot_silhouette()}, \code{model$summary()}
#'   \item Predict new observations: \code{predictions <- model$predict(new_data)}
#' }
#'
#' @section Dependencies:
#' Required packages:
#' \itemize{
#'   \item \pkg{R6}: For the R6 object-oriented class system
#'   \item \pkg{FactoMineR}: For creating disjunctive tables and performing MCA
#' }
#' Optional packages:
#' \itemize{
#'   \item \pkg{dendextend}: For enhanced dendrogram visualization with colored branches
#' }
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(R6)
#' library(FactoMineR)
#'
#' # Create sample categorical data
#' set.seed(123)
#' data <- data.frame(
#'   color = factor(sample(c("red", "blue", "green"), 100, replace = TRUE)),
#'   size = factor(sample(c("small", "medium", "large"), 100, replace = TRUE)),
#'   shape = factor(sample(c("circle", "square", "triangle"), 100, replace = TRUE))
#' )
#'
#' # Approach 1: Specify number of clusters at initialization
#' model <- ModCluster$new(n_clusters = 3, hclust_method = "average")
#' model$fit(data)
#'
#' # Approach 2: Fit first, then decide on number of clusters
#' model <- ModCluster$new()
#' model$fit(data)
#' model$plot_heights()  # Examine elbow plot
#' model$cut_tree(k = 3)  # Cut tree at desired level
#'
#' # Visualize clustering results
#' model$plot_dendrogram()  # Hierarchical tree
#' model$plot_mca()  # MCA factorial plane
#' model$plot_silhouette()  # Clustering quality
#'
#' # Get clustering results
#' model$print()  # Concise summary
#' cluster_table <- model$get_cluster_table()  # Simple table
#' summary_results <- model$summary()  # Detailed statistics
#'
#' # Access clustering information
#' print(summary_results$clust_summary)  # Cluster sizes
#' print(summary_results$clust_members)  # Modality details
#'
#' # Create new data for prediction
#' new_data <- data.frame(
#'   color = factor(c("red", "blue"), levels = c("red", "blue", "green")),
#'   size = factor(c("small", "large"), levels = c("small", "medium", "large")),
#'   shape = factor(c("circle", "square"), levels = c("circle", "square", "triangle"))
#' )
#'
#' # Predict cluster membership
#' predictions <- model$predict(new_data)
#' print(predictions)
#'
#' # Visualize new observations in MCA space
#' model$plot_mca_with_new(new_data)
#'
#' # Project new observations onto MCA axes
#' projected <- model$project_new_modalities(new_data, axes = c(1, 2))
#' print(projected)
#'
#' # Access active bindings
#' labels <- model$modality_labels  # Cluster assignments
#' freqs <- model$modality_frequencies  # Modality counts
#' centers <- model$cluster_centers  # Cluster profiles
#' }
#'
#' @seealso
#' \code{\link[FactoMineR]{MCA}} for Multiple Correspondence Analysis,
#' \code{\link[stats]{hclust}} for hierarchical clustering,
#' \code{\link[stats]{dist}} for distance computation
#'
#' @references
#' Dice, L. R. (1945). Measures of the amount of ecologic association between species.
#' \emph{Ecology}, 26(3), 297-302.
#'
#' @name ModCluster
#' @aliases modality_clustering dice_clustering categorical_clustering
#' @export
ModCluster <- R6::R6Class(
  "ModCluster",
  inherit = ClusteringBase,

  private = list(
    .method = "hclust",
    .disjunctive_matrix = NULL,
    .disjunctive_matrix_illust = NULL,
    .data_illust = NULL,
    .dice_matrix = NULL,
    .modality_dist = NULL,
    .clustering_object = NULL,
    .n_dimensions = NULL,
    .modality_names = NULL,
    .modality_labels = NULL,
    .hclust_method = NULL,
    .modality_frequencies = NULL,
    .summary_results = NULL,
    .mca_result = NULL,
    .modality_coords = NULL,
    .variable_info = NULL,

    #' Build the complete disjunctive matrix with unique modality names
    #' @description Internal method to create indicator matrix from categorical data
    #' @return Matrix with binary indicators for each modality
    create_disjunctive_matrix = function() {
      if (length(private$.quali_indices) == 0) {
        stop("No qualitative variables found. This algorithm requires categorical data.")
      }

      data_quali <- self$get_quali_data()
      disj_matrix <- FactoMineR::tab.disjonctif(data_quali)

      # ALWAYS create unique names in format variable.level
      var_names <- names(data_quali)
      new_names <- character(ncol(disj_matrix))
      col_idx <- 1

      for (var_name in var_names) {
        levels_var <- levels(data_quali[[var_name]])
        n_levels <- length(levels_var)
        new_names[col_idx:(col_idx + n_levels - 1)] <-
          paste(var_name, levels_var, sep = ".")
        col_idx <- col_idx + n_levels
      }

      colnames(disj_matrix) <- new_names
      private$.modality_names <- new_names
      private$.modality_frequencies <- colSums(disj_matrix)

      # Save variable info for predict()
      private$.variable_info <- list(
        var_names = var_names,
        var_levels = lapply(data_quali, levels)
      )

      message("✓ Disjunctive matrix created: ",
              nrow(disj_matrix), " individuals × ",
              ncol(disj_matrix), " modalities")

      return(disj_matrix)
    },

    #' Compute Dice distances between all modalities
    #' @description Internal method to calculate pairwise Dice distances
    #' @return NULL (stores result in private$.dice_matrix)
    compute_dice_distances = function() {
      disj <- private$.disjunctive_matrix

      n <- nrow(disj)
      m <- ncol(disj)

      dice_matrix <- matrix(0, nrow = m, ncol = m)
      message("  Computing Dice distances...")

      for (j in 1:m) {
        for (j_prime in j:m) {
          if (j == j_prime) {
            dice_matrix[j, j_prime] <- 0
          } else {
            val <- 0.5 * sum((disj[, j] - disj[, j_prime])**2)

            dice_matrix[j, j_prime] <- val
            dice_matrix[j_prime, j] <- val
          }
        }
      }

      colnames(dice_matrix) <- private$.modality_names
      rownames(dice_matrix) <- private$.modality_names
      private$.dice_matrix <- dice_matrix
      private$.modality_dist <- as.dist(dice_matrix)

      message("✓ Dice distance matrix computed")
      message("  - ", m, " modalities")
      message("  - ", m*(m-1)/2, " pairwise distances")

      return(invisible(NULL))
    },

    #' Compute MCA for visualization
    #' @description Internal method to perform Multiple Correspondence Analysis
    #' @return NULL (stores result in private$.mca_result)
    compute_mca_for_visualization = function() {
      data_quali <- self$get_quali_data()

      max_possible <- ncol(private$.disjunctive_matrix) - ncol(data_quali)
      n_dims <- min(private$.n_dimensions, max_possible)

      mca_res <- FactoMineR::MCA(
        data_quali,
        ncp = n_dims,
        graph = FALSE
      )

      private$.mca_result <- mca_res
      private$.modality_coords <- mca_res$var$coord

      # Ensure MCA row names match disjunctive matrix column names
      var_names <- names(data_quali)
      levels_list <- lapply(data_quali, levels)

      mca_names <- unlist(mapply(
        function(v, lv) paste(v, lv, sep = "."),
        var_names,
        levels_list,
        SIMPLIFY = FALSE
      ))

      if (length(mca_names) != nrow(private$.modality_coords)) {
        warning("Mismatch between MCA coordinates and expected modality names")
      }

      rownames(private$.mca_result$var$coord) <- mca_names
      rownames(private$.modality_coords) <- mca_names

      if (n_dims < private$.n_dimensions) {
        message("  Note: MCA dimensions reduced to ", n_dims)
        private$.n_dimensions <- n_dims
      }

      message("✓ MCA computed (visualization only)")
      message("  - ", n_dims, " dimensions")
      message("  - Variance explained: ",
              round(sum(mca_res$eig[1:n_dims, 2]), 2), "%")

      return(invisible(NULL))
    },

    #' Apply hierarchical clustering
    #' @description Internal method to perform hierarchical clustering on Dice distances
    #' @return NULL (stores result in private$.clustering_object)
    apply_clustering = function() {
      if (is.null(private$.modality_dist)) {
        stop("Dice distances must be computed first")
      }

      hc <- hclust(private$.modality_dist, method = private$.hclust_method)
      private$.clustering_object <- hc

      if (!is.null(private$.n_clusters)) {
        private$.modality_labels <- cutree(hc, k = private$.n_clusters)
        names(private$.modality_labels) <- private$.modality_names

        message("✓ Hierarchical clustering completed")
        message("  - ", private$.n_clusters, " clusters")
        message("  - Linkage: ", private$.hclust_method)
      } else {
        message("✓ Hierarchical clustering completed")
        message("  - Use cut_tree(k) to specify number of clusters")
      }
    },



    #' Create disjunctive matrix for new data
    #' @description Internal method to create indicator matrix for prediction
    #' @param new_data Data frame with same variables as training data
    #' @return Matrix with binary indicators
    create_new_disjunctive = function(new_data) {
      if (is.null(private$.variable_info)) {
        stop("Model must be fitted first")
      }

      var_info <- private$.variable_info
      var_names <- var_info$var_names

      missing_vars <- setdiff(var_names, names(new_data))
      if (length(missing_vars) > 0) {
        stop("Missing variables in new_data: ", paste(missing_vars, collapse = ", "))
      }

      new_data_quali <- new_data[, var_names, drop = FALSE]
      for (i in seq_along(var_names)) {
        var <- var_names[i]
        expected_levels <- var_info$var_levels[[i]]
        new_data_quali[[var]] <- factor(new_data_quali[[var]], levels = expected_levels)

        unknown_levels <- setdiff(unique(new_data[[var]]), expected_levels)
        if (length(unknown_levels) > 0) {
          warning("Unknown levels in variable '", var, "': ",
                  paste(unknown_levels, collapse = ", "))
        }
      }

      disj_new <- FactoMineR::tab.disjonctif(new_data_quali)

      expected_names <- private$.modality_names
      if (!all(colnames(disj_new) %in% expected_names)) {
        disj_aligned <- matrix(0, nrow = nrow(disj_new), ncol = length(expected_names))
        colnames(disj_aligned) <- expected_names
        common_cols <- intersect(colnames(disj_new), expected_names)
        disj_aligned[, common_cols] <- disj_new[, common_cols]
        disj_new <- disj_aligned
      }

      return(disj_new)
    },

    #' Compute Dice distance to modalities
    #' @description Internal method to calculate distances from new profile to existing modalities
    #' @param new_profile Binary vector representing a new observation
    #' @return Named vector of distances
    compute_dice_to_modalities = function(new_profile) {
      disj_train <- private$.disjunctive_matrix
      n_train <- nrow(disj_train)
      m <- ncol(disj_train)

      dice_dists <- numeric(m)

      for (j in 1:m) {
        val <- 0.5 * sum(
          new_profile * (1 - disj_train[, j]) +
            (1 - new_profile) * disj_train[, j]
        )
        dice_dists[j] <- sqrt(val)
      }

      names(dice_dists) <- private$.modality_names
      return(dice_dists)
    },

    #' @description
    #' Project new observations into MCA space
    #' @param axes Integer vector of MCA axes to use (default: all)
    #' @return Matrix of projected coordinates
    project_new_modalities = function(axes = NULL) {

      # Check if fitted
      if (!private$.fitted) {
        stop("Model must be fitted")
      }

      # Get the illustrative data
      new_data <- private$.data_illust

      # Check if predict have been made.
      if (is.null(new_data)) {
        stop("Make a prediction first")
      }

      message("Projecting new modalities into MCA space...")

      # Get the MCA result
      mca_fit <- private$.mca_result

      # Get individual coordinates from FactoMineR MCA result
      ind_coords <- mca_fit$ind$coord

      # Get eigenvalues
      eigenvalues <- mca_fit$eig[, 1]  # First column contains eigenvalues


      # Determine axes to use
      if (is.null(axes)) {
        axes <- 1:ncol(ind_coords)
      }

      # For each new variable: compute mean coordinates for each modality
      new_mod_coords_list <- list()
      n_new <- 0

      for (var_name in colnames(new_data)) {
        # Get the variable values
        var_values <- new_data[[var_name]]

        # For each unique modality of this variable
        unique_mods <- unique(var_values)
        unique_mods <- unique_mods[!is.na(unique_mods)]  # Remove NAs

        for (mod in unique_mods) {
          # Identify individuals with this modality
          idx <- which(var_values == mod)
          n_new <- n_new + length(idx)

          # Compute mean of factorial coordinates for these individuals
          mean_coords <- colMeans(ind_coords[idx, axes, drop = FALSE])

          # Divide by square root of eigenvalues to get true coordinates
          mean_coords <- mean_coords / sqrt(eigenvalues[axes])

          # Store with label "variable.modality"
          mod_label <- paste0(var_name, ".", mod)
          new_mod_coords_list[[mod_label]] <- mean_coords
        }
      }

      # Convert to matrix
      new_mod_coords <- do.call(rbind, new_mod_coords_list)
      colnames(new_mod_coords) <- paste0("Dim ", axes)

      message("Projection completed for ", nrow(new_data), " observations")
      message("  New modalities projected: ", length(new_mod_coords_list))
      message("  Axes used: ", paste(axes, collapse = ", "))

      return(new_mod_coords)
    }
  ),

  public = list(
    #' @description
    #' Initialize a new ModCluster object
    #' @param n_clusters Integer, target number of clusters (default: NULL)
    #' @param method Character, clustering method (default: "hclust")
    #' @param n_dimensions Integer, number of MCA dimensions (default: 5)
    #' @param hclust_method Character, linkage method (default: "average")
    #' @return A new ModCluster object
    initialize = function(method = "hclust",
                          n_dimensions = 2,
                          hclust_method = "average") {

      if (method != "hclust") {
        stop("Only 'hclust' is supported with Dice distances")
      }

      valid_methods <- c("single", "complete", "average")

      if (!hclust_method %in% valid_methods) {
        stop("Invalid hclust_method. Choose from: ",
             paste(valid_methods, collapse = ", "))
      }

      private$.method <- method
      private$.n_dimensions <- n_dimensions
      private$.hclust_method <- hclust_method

      message("  Linkage: ", hclust_method)
      message("  MCA dims (viz): ", n_dimensions)
    },


    #' @description
    #' Fit the clustering model on training data
    #' @param data Data frame with categorical variables
    #' @return Self (invisibly) for method chaining
    fit = function(data) {
      self$load_and_check_data(data)
      self$validate_algorithm_requirements("qual")

      message("\n[1/4] Building disjunctive matrix...")
      private$.disjunctive_matrix <- private$create_disjunctive_matrix()

      message("\n[2/4] Computing Dice distances...")
      private$compute_dice_distances()

      message("\n[3/4] Applying hierarchical clustering...")
      private$apply_clustering()

      message("\n[4/4] Computing MCA for visualization...")
      private$compute_mca_for_visualization()

      self$fitted <- TRUE
      message("\n✔ Modality clustering completed successfully!\n")

      return(invisible(self))
    },

    #' @description
    #' Predict cluster membership for new observations
    #' @param new_data Data frame with same variables as training data
    #' @return Data frame with predictions and distances
    predict = function(new_data) {

      if (!self$fitted) {
        stop("Model must be fitted")
      }

      if (nrow(new_data) != nrow(private$.data)) {
        stop("New values have invalid row numbers.")

      } else if (is.null(private$.modality_labels)) {
        stop("Tree must be cutted first")
      }

      quant_cols <- sapply(new_data, is.numeric)
      if (any(quant_cols)) {
        warning("Quantitative columns discarded : ",
                paste(names(new_data)[quant_cols], collapse = ", "))
      }

      # Keep only non quantitatives
      new_data <- new_data[ , !quant_cols, drop = FALSE]

      # Convert to factor
      new_data[] <- lapply(new_data, as.factor)

      # Add the illustrative data into the private field
      private$.data_illust <- new_data

      # transform into disjunctive table.
      disj_matrix_illust <- FactoMineR::tab.disjonctif(new_data)


      # Create new modalities name considering the variable.
      var_names <- colnames(new_data)
      new_names <- character(ncol(disj_matrix_illust))
      col_idx <- 1

      data_quali <- self$get_quali_data()

      # Set the right name to the new ilustrative variables. (format newVariable.modality).
      for (var_name in var_names) {

        # Get the moldalities of the given illust variable
        levels_var <- levels(new_data[[var_name]])

        # Get the number of modalities
        n_levels <- length(levels_var)

        # Copy the modality vector into the new_names vector
        new_names[col_idx:(col_idx + n_levels - 1)] <-
          paste(var_name, levels_var, sep = ".")
        col_idx <- col_idx + n_levels
      }

      colnames(disj_matrix_illust) <- new_names

      # Add the disjunctive matrix illustrative to private field
      private$.disjunctive_matrix_illust <- disj_matrix_illust

      # Calculate the dice matrix of illustrative variable with all data.
      d2 <- matrix(0, ncol(disj_matrix_illust), ncol(private$.disjunctive_matrix))
      for (j in 1:ncol(disj_matrix_illust)) {
        for(jprim in 1:ncol(private$.disjunctive_matrix)) {
          d2[j, jprim] <- 0.5 * (sum((disj_matrix_illust[, j] - private$.disjunctive_matrix[, jprim])**2))
        }
      }
      colnames(d2) <- colnames(private$.disjunctive_matrix)
      rownames(d2) <- colnames(disj_matrix_illust)
      # Initialize cluster_distances matrix.
      cluster_distances <- data.frame(
        modality = colnames(disj_matrix_illust)
      )
      # Calculate cluster distances.
      for (i in 1:private$.n_clusters) {
        # Get indices of modalities in cluster i
        cluster_indices <- which(private$.modality_labels == i)

        # Initialize the vector of distances of each variable to the ith cluster.
        temp_distances <- c()

        if (private$.hclust_method == "average") {

          # Compute the mean distance between each illustrative variable and a cluster
          if (length(cluster_indices) == 1) {

            # If only one modality in cluster, use it directly
            temp_distances <- d2[, cluster_indices]
          } else {
            temp_distances <- apply(d2[, cluster_indices, drop = FALSE], MARGIN = 1, FUN = mean)
          }
        } else if (private$.hclust_method == "single"){

          # compute the minimum distance
          if (length(cluster_indices) == 1) {

            temp_distances <- d2[, cluster_indices]
          } else {
            temp_distances <- apply(d2[, cluster_indices, drop = FALSE], MARGIN = 1, FUN = min)
          }
        } else if (private$.hclust_method == "complete") {

          # compute the maximum distance to clusters
          if (length(cluster_indices) == 1) {
            temp_distances <- d2[, cluster_indices]
          } else {

            temp_distances <- apply(d2[, cluster_indices, drop = FALSE], MARGIN = 1, FUN = max)
          }
        }
        # Save the temp distances vector to the cluster_distances dataframe
        cluster_distances[[paste0("cluster_", i)]] <- temp_distances
      }

      # Cluster attribution by the minimal distance.
      predicted_labels <- apply(cluster_distances[, -1], 1, which.min)
      names(predicted_labels) <- colnames(disj_matrix_illust)

      # return a list of object containing the new labels and distances
      return(list(prediction = predicted_labels, distances = cluster_distances))
    },


    #' @description
    #' Get the Dice distance matrix
    #' @return Symmetric matrix of Dice distances
    get_dice_matrix = function() {
      if (!self$fitted) stop("Model must be fitted first")
      return(private$.dice_matrix)
    },

    #' @description
    #' Print a summary of the ModCluster object
    #' @return Self (invisibly)
    print = function() {
      if (!self$fitted) {
        cat("Status: NOT FITTED\n")
        cat("Use $fit(data) to train the model.\n\n")
        return(invisible(self))
      }

      cat("Method: Hierarchical Clustering (", private$.hclust_method, ")\n", sep = "")
      cat("Modalities: ", length(private$.modality_names), "\n", sep = "")

      if (!is.null(private$.n_clusters)) {
        cat("Clusters: ", private$.n_clusters, "\n", sep = "")
        cat("\nCluster sizes:\n")
        print(table(private$.modality_labels))
      } else {
        cat("Clusters: Not yet specified (use cut_tree)\n")
      }

      invisible(self)
    },

    #' @description
    #' Get a table of modalities with cluster assignments
    #' @return Data frame with Modality, Cluster, and Frequency columns
    get_cluster_table = function() {
      if (!self$fitted) stop("Model must be fitted")
      if (is.null(private$.modality_labels)) {
        stop("No clusters defined. Use cut_tree(k) first.")
      }

      df <- data.frame(
        Modality = names(private$.modality_labels),
        Cluster = private$.modality_labels,
        Frequency = private$.modality_frequencies[names(private$.modality_labels)],
        stringsAsFactors = FALSE
      )

      df <- df[order(df$Cluster, -df$Frequency), ]
      rownames(df) <- NULL
      return(df)
    },

    #' @description
    #' Cut the hierarchical tree into k clusters
    #' @param k Integer, number of clusters
    #' @return Self (invisibly) for method chaining
    cut_tree = function(k) {
      if (!self$fitted) stop("Model must be fitted before cut_tree()")

      private$.n_clusters <- k
      private$.modality_labels <- cutree(private$.clustering_object, k = k)
      names(private$.modality_labels) <- private$.modality_names
      private$.summary_results <- NULL

      message("✓ Tree cut into ", k, " clusters")
      invisible(self)
    },

    #' @description
    #' Plot hierarchical clustering dendrogram
    #' @param k Integer, number of clusters to highlight (default: current k)
    #' @param color_branches Logical, color branches by cluster (default: TRUE)
    #' @return NULL (invisibly)
    plot_dendrogram = function(k = NULL, color_branches = TRUE) {
      if (!self$fitted) stop("Model must be fitted first")

      if (is.null(k)) k <- private$.n_clusters

      if (is.null(k)) {
        plot(private$.clustering_object,
             main = "Dendrogram - Modality Clustering (Dice Distance)",
             xlab = "Modalities",
             ylab = "Height",
             sub = paste("Linkage:", private$.hclust_method))
        return(invisible(NULL))
      }

      if (color_branches && requireNamespace("dendextend", quietly = TRUE)) {
        dend <- as.dendrogram(private$.clustering_object)
        dend_colored <- dendextend::color_branches(dend, k = k)
        plot(dend_colored,
             main = paste("Dendrogram -", k, "Clusters (Dice Distance)"),
             xlab = "Modalities",
             ylab = "Height")
      } else {
        plot(private$.clustering_object,
             main = paste("Dendrogram -", k, "Clusters (Dice Distance)"),
             xlab = "Modalities",
             ylab = "Height")
        rect.hclust(private$.clustering_object, k = k, border = 2:6)
      }
    },

    #' @description
    #' Plot modalities in MCA space
    #' @param axes Integer vector of length 2, which axes to plot (default: c(1, 2))
    #' @param show_labels Logical, show modality labels (default: TRUE)
    #' @param label_cex Numeric, label text size (default: 0.7)
    #' @return NULL (invisibly)
    plot_mca = function(axes = c(1, 2), show_labels = TRUE, label_cex = 0.7) {
      if (!self$fitted) stop("Model must be fitted first")
      if (is.null(private$.mca_result)) {
        message("Computing MCA...")
        private$compute_mca_for_visualization()
      }
      if (is.null(private$.modality_labels)) {
        stop("No clusters defined. Use cut_tree(k) first.")
      }
      coords <- private$.modality_coords
      labels <- private$.modality_labels
      common_names <- intersect(rownames(coords), names(labels))
      if (length(common_names) == 0) {
        cat("ERROR - MCA modalities:\n")
        print(rownames(coords))
        cat("\nERROR - Clustering modalities:\n")
        print(names(labels))
        stop("No common modalities between MCA and clustering results")
      }
      coords <- coords[common_names, , drop = FALSE]
      labels <- labels[common_names]
      coords2d <- coords[, axes, drop = FALSE]

      # Calculate symmetric limits around 0
      max_x <- max(abs(coords2d[, 1])) * 1.1
      max_y <- max(abs(coords2d[, 2])) * 1.1
      xlim <- c(-max_x, max_x)
      ylim <- c(-max_y, max_y)

      colors <- rainbow(private$.n_clusters, alpha = 0.8)[labels]
      eig <- private$.mca_result$eig
      xlab <- sprintf("Dim %d (%.1f%%)", axes[1], eig[axes[1], 2])
      ylab <- sprintf("Dim %d (%.1f%%)", axes[2], eig[axes[2], 2])

      plot(
        coords2d,
        col = colors,
        pch = 19,
        cex = 1.2,
        xlab = xlab,
        ylab = ylab,
        main = paste("MCA Projection -", private$.n_clusters, "Clusters"),
        xlim = xlim,
        ylim = ylim,
        panel.first = {
          grid(col = "gray90", lty = 1)
          abline(h = 0, v = 0, lty = 2, col = "gray50")
        }
      )
      if (show_labels) {
        text(coords2d, labels = rownames(coords2d),
             pos = 3, cex = label_cex, col = colors)
      }
      legend(
        "topright",
        legend = paste("Cluster", 1:private$.n_clusters),
        col = rainbow(private$.n_clusters, alpha = 0.8),
        pch = 19,
        cex = 0.8,
        bg = "white"
      )
      invisible(NULL)
    },

    #' @description
    #' Plot training modalities (points) and new observations in MCA space
    #' @param axes Integer vector of length 2, which axes to plot (default: c(1, 2))
    #' @param show_labels Logical, show labels (default: TRUE)
    #' @param label_cex Numeric, label text size (default: 0.7)
    #' @return Invisible list with coordinates and cluster labels
    plot_mca_with_new = function(axes = c(1, 2), show_labels = TRUE, label_cex = 0.7) {
      # Vérifications préliminaires
      if (!self$fitted) stop("Model must be fitted first")
      if (is.null(private$.mca_result)) {
        message("Computing MCA...")
        private$compute_mca_for_visualization()
      }
      if (is.null(private$.modality_labels)) {
        stop("No clusters defined. Use cut_tree(k) first.")
      }
      # 1. Get training modality coordinates
      coords_train <- private$.modality_coords[, axes, drop = FALSE]
      labels_train <- private$.modality_labels
      eig <- private$.mca_result$eig
      # 2. Project new modalities
      coords_new <- private$project_new_modalities(axes = axes)
      # 3. Extract cluster labels for new modalities
      labels_new <- as.numeric(gsub(".*_cluster", "", rownames(coords_new)))
      # 4. Prepare colors
      color_palette <- rainbow(private$.n_clusters, alpha = 0.8)
      colors_train <- color_palette[labels_train]
      colors_new <- color_palette[labels_new]
      # 5. Calculate plot limits (SYMMETRIC around 0)
      all_coords <- rbind(coords_train, coords_new)
      max_x <- max(abs(all_coords[, 1])) * 1.1
      max_y <- max(abs(all_coords[, 2])) * 1.1
      xlim <- c(-max_x, max_x)
      ylim <- c(-max_y, max_y)
      # 6. Create plot
      variance_explained <- eig[axes, 2]
      plot(coords_train[, 1], coords_train[, 2],
           xlim = xlim, ylim = ylim,
           xlab = sprintf("Dim %d (%.2f%%)", axes[1], variance_explained[1]),
           ylab = sprintf("Dim %d (%.2f%%)", axes[2], variance_explained[2]),
           main = "MCA Factor Map: Training vs New Modalities",
           col = colors_train,
           pch = 16,
           cex = 1.2)
      # Add new modalities with different symbol
      points(coords_new[, 1], coords_new[, 2],
             col = colors_new,
             pch = 17,  # triangles for new modalities
             cex = 1.5)
      # Add axes (now centered)
      abline(h = 0, v = 0, col = "gray", lty = 2, lwd = 1.5)
      # Add labels if requested
      if (show_labels) {
        text(coords_train[, 1], coords_train[, 2],
             labels = rownames(coords_train),
             col = colors_train,
             cex = label_cex,
             pos = 3,
             offset = 0.3)
        text(coords_new[, 1], coords_new[, 2],
             labels = rownames(coords_new),
             col = colors_new,
             cex = label_cex,
             pos = 1,
             offset = 0.3,
             font = 2)  # bold for new modalities
      }
      # Add legend
      legend("topright",
             legend = c(paste("Cluster", 1:private$.n_clusters),
                        "Training", "New"),
             col = c(color_palette, "black", "black"),
             pch = c(rep(16, private$.n_clusters), 16, 17),
             bty = "n",
             cex = 0.8)
      # Return coordinates invisibly
      result <- list(
        coords_train = coords_train,
        coords_new = coords_new,
        variance = variance_explained
      )
      invisible(result)
    },

    #' @description
    #' Plot silhouette values for modality clustering
    #' @return Invisible data frame with silhouette values
    plot_silhouette = function() {
      if (!self$fitted) stop("Model must be fitted first")
      if (is.null(private$.modality_labels)) {
        stop("No clusters defined. Use cut_tree(k) first.")
      }

      n_mod <- length(private$.modality_labels)
      labels <- private$.modality_labels
      dice_mat <- private$.dice_matrix

      sil_values <- numeric(n_mod)

      for (i in 1:n_mod) {
        mod_cluster <- labels[i]
        same_cluster <- which(labels == mod_cluster & names(labels) != names(labels)[i])

        if (length(same_cluster) > 0) {
          a_i <- mean(dice_mat[i, same_cluster])
        } else {
          a_i <- 0
        }

        other_clusters <- setdiff(unique(labels), mod_cluster)
        min_dist_to_other <- Inf

        for (cl in other_clusters) {
          other_members <- which(labels == cl)
          if (length(other_members) > 0) {
            d <- mean(dice_mat[i, other_members])
            min_dist_to_other <- min(min_dist_to_other, d)
          }
        }
        b_i <- min_dist_to_other

        if (b_i > 0) {
          sil <- (b_i-a_i)/max(a_i,b_i)
          sil_values[i] <- sil
        } else {
          sil_values[i] <- 0
        }
      }

      sil_df <- data.frame(
        modality = names(labels),
        cluster = labels,
        silhouette = sil_values,
        stringsAsFactors = FALSE
      )

      sil_df <- sil_df[order(sil_df$cluster, sil_df$silhouette), ]

      col_palette <- rainbow(private$.n_clusters)
      cluster_colors <- col_palette[sil_df$cluster]

      y_pos <- 1:nrow(sil_df)
      x_range <- range(c(0, sil_df$silhouette))

      par(mar = c(5, 8, 4, 2))
      plot(
        sil_df$silhouette,
        y_pos,
        type = "n",
        xlab = "Silhouette value ",
        ylab = "",
        yaxt = "n",
        main = "Silhouette Plot - Modality Clustering",
        xlim = x_range,
        ylim = c(0.5, nrow(sil_df) + 0.5)
      )

      abline(v = 0, col = "gray70", lty = 2)
      abline(v = mean(sil_df$silhouette), col = "red", lty = 2, lwd = 2)

      for (i in 1:nrow(sil_df)) {
        rect(
          xleft = 0,
          xright = sil_df$silhouette[i],
          ybottom = y_pos[i] - 0.4,
          ytop = y_pos[i] + 0.4,
          col = cluster_colors[i],
          border = NA
        )
      }

      cluster_changes <- which(diff(sil_df$cluster) != 0)
      if (length(cluster_changes) > 0) {
        abline(h = y_pos[cluster_changes] + 0.5, col = "black", lwd = 2)
      }

      cluster_mids <- tapply(y_pos, sil_df$cluster, mean)
      axis(2, at = cluster_mids, labels = paste("C", names(cluster_mids)), las = 1)

      avg_sil <- mean(sil_df$silhouette)
      mtext(sprintf("Average silhouette: %.3f", avg_sil), side = 3, line = 0.5)

      legend(
        "topright",
        legend = c(paste("Cluster", 1:private$.n_clusters), "Mean"),
        fill = c(col_palette, NA),
        border = c(rep("black", private$.n_clusters), NA),
        lty = c(rep(NA, private$.n_clusters), 2),
        lwd = c(rep(NA, private$.n_clusters), 2),
        col = c(rep(NA, private$.n_clusters), "red"),
        cex = 0.7,
        bg = "white"
      )

      par(mar = c(5, 4, 4, 2))

      invisible(sil_df)
    },

    #' @description
    #' Plot aggregation heights (elbow plot)
    #' @return NULL (invisibly)
    plot_heights = function() {
      if (!self$fitted) {
        stop("Model must be fitted first")
      }

      heights <- private$.clustering_object$height
      n <- length(heights) + 1
      k_values <- 2:(n - 1)
      heights_k <- heights[n - k_values]


      plot(
        k_values,
        heights_k,
        type = "b",
        pch = 19,
        xlab = "Number of clusters",
        ylab = "Aggregation height",
        main = "Aggregation Heights (Elbow Plot)",
        panel.first = grid(col = "gray90")
      )

      if (!is.null(private$.n_clusters)) {
        k_current <- private$.n_clusters
        if (k_current >= 2 && k_current <= (n - 1)) {
          abline(v = k_current, col = "red", lwd = 2, lty = 2)
          text(k_current, max(heights_k) * 0.9,
               labels = paste("k =", k_current),
               col = "red", pos = 4)
        }
      }

      invisible(NULL)
    },

    #' @description
    #' Get comprehensive clustering statistics
    #' @return List with cluster summary and member details
    summary = function() {
      if (!self$fitted || is.null(private$.modality_labels)) {
        cat("No clusters defined. Use cut_tree(k) first.\n")
        return(invisible(NULL))
      }

      if (!is.null(private$.summary_results)) {
        return(private$.summary_results)
      }

      n_clust <- private$.n_clusters
      labels <- private$.modality_labels
      dice_mat <- private$.dice_matrix

      clust_summary <- data.frame(
        cluster = 1:n_clust,
        n_members = as.vector(table(labels))
      )

      n_mod <- length(labels)
      clust_members <- data.frame(
        modality = names(labels),
        cluster = labels,
        frequency = private$.modality_frequencies[names(labels)],
        stringsAsFactors = FALSE
      )

      own_dist <- numeric(n_mod)
      next_dist <- numeric(n_mod)

      for (i in 1:n_mod) {
        mod_cluster <- labels[i]
        same_cluster <- which(labels == mod_cluster & names(labels) != names(labels)[i])

        if (length(same_cluster) > 0) {
          if (private$.hclust_method == "average") {
            own_dist[i] <- mean(dice_mat[i, same_cluster])
          } else if (private$.hclust_method == "single") {
            own_dist[i] <- min(dice_mat[i, same_cluster])
          } else if (private$.hclust_method == "complete") {
            own_dist[i] <- max(dice_mat[i, same_cluster])
          } else {
            own_dist[i] <- mean(dice_mat[i, same_cluster])
          }
        } else {
          own_dist[i] <- 0
        }

        other_clusters <- setdiff(1:n_clust, mod_cluster)
        min_dist_to_other <- Inf

        for (cl in other_clusters) {
          other_members <- which(labels == cl)
          if (length(other_members) > 0) {
            if (private$.hclust_method == "average") {
              d <- mean(dice_mat[i, other_members])
            } else if (private$.hclust_method == "single") {
              d <- min(dice_mat[i, other_members])
            } else {
              d <- max(dice_mat[i, other_members])
            }
            min_dist_to_other <- min(min_dist_to_other, d)
          }
        }
        next_dist[i] <- min_dist_to_other
      }

      clust_members$dist_to_own <- round(sqrt(own_dist), 3)
      clust_members$dist_to_next <- round(sqrt(next_dist), 3)
      clust_members$ratio <- round(clust_members$dist_to_own / clust_members$dist_to_next, 3)

      clust_members <- clust_members[order(clust_members$cluster, clust_members$dist_to_own), ]
      rownames(clust_members) <- NULL

      result <- list(
        clust_summary = clust_summary,
        clust_members = clust_members
      )

      private$.summary_results <- result
      return(result)
    }
  ),

  active = list(
    #' @field modality_labels Named integer vector of cluster assignments
    modality_labels = function() private$.modality_labels,

    #' @field modality_names Character vector of modality names
    modality_names = function() private$.modality_names,

    #' @field modality_frequencies Named numeric vector of modality frequencies
    modality_frequencies = function() private$.modality_frequencies,

    #' @field mca_result Complete MCA result object
    mca_result = function() private$.mca_result,

    #' @field data_illust Return the illustrative data.
    data_illust = function() private$.data_illust
  )
)

#' @title Print Method for ModCluster
#' @description Print a summary of the ModCluster object
#' @param x A ModCluster object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the object
#' @export
print.ModCluster <- function(x, ...) {
  x$print()
  invisible(x)
}

#' @title Summary Method for ModCluster
#' @description Get detailed clustering statistics
#' @param object A ModCluster object
#' @param ... Additional arguments (ignored)
#' @return A list containing cluster summary and member details
#' @export
summary.ModCluster <- function(object, ...) {
  object$summary()
}
