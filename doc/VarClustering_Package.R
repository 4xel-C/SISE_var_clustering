## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 10,
  fig.height = 6,
  fig.align = "center"
)

## ----libraries----------------------------------------------------------------
library(VarClustering)

## ----hclust-load-auto---------------------------------------------------------
# Load the autos dataset
data(autos)

# Display structure
str(autos)

# Summary
summary(autos)

## ----hclust-split-auto--------------------------------------------------------
# Split into active and illustrative variables
autos_train <- autos[, 2:9]
autos_test <- autos[, 10:12]

## ----hclust-fit-auto----------------------------------------------------------
# Initialize HClustVar
hc_auto <- HClustVar$new(
  vartype = "quant",
  dist.metric = "rsquare",
  cah.method = "ward.D2"
)

# Fit the model
hc_auto$fit(autos_train)

# Print summary
print(hc_auto)

## ----hclust-dendro-auto-------------------------------------------------------
hc_auto$plot_dendrogram()

## ----hclust-elbow-auto--------------------------------------------------------
hc_auto$plot_agg_levels()

## ----hclust-cut-auto----------------------------------------------------------
# Cut into 3 clusters
hc_auto$cut_tree(k = 3)

# Visualize with cluster colors
hc_auto$plot_dendrogram(k = 3)

## ----hclust-quality-auto------------------------------------------------------
# MDS projection
hc_auto$mds_projection()

# Silhouette plot
hc_auto$plot_silhouette()

## ----hclust-summary-auto------------------------------------------------------
results <- hc_auto$summary()

# Quality metrics
cat("Average silhouette:", results$avg_silhouette, "\n")
cat("Variance explained:", results$total_var_explained, "\n")

# Cluster summary
results$clust_summary

## ----hclust-members-auto------------------------------------------------------
# Cluster members
results$clust_members

## ----hclust-corr-auto---------------------------------------------------------
# Centroid correlations
results$centroids_correlations

## ----hclust-predict-auto------------------------------------------------------
# Predict cluster for illustrative variables
prediction <- hc_auto$predict(autos_test)

# Labels
prediction$labels

# Proximities to all clusters
prediction$proximities

## ----hclust-canines-----------------------------------------------------------
# Load canines dataset
data(canines)

# Display
summary(canines)
head(canines)

## ----hclust-fit-canines-------------------------------------------------------
# Select active variables
canines_actives <- canines[, 2:7]

# Initialize for qualitative variables
hc_canines <- HClustVar$new(
  vartype = "qual",
  cah.method = "ward.D2"
)

# Fit
hc_canines$fit(canines_actives)

# Dendrogram
hc_canines$plot_dendrogram()

## ----hclust-elbow-canines-----------------------------------------------------
# Aggregation levels
hc_canines$plot_agg_levels()

## ----hclust-cut-canines-------------------------------------------------------
# Cut and analyze
hc_canines$cut_tree(4)
canines_result <- hc_canines$summary()

# Results
canines_result$clust_summary
canines_result$clust_members

## ----hclust-mds-canines-------------------------------------------------------
# MDS visualization
hc_canines$mds_projection()

## ----kmeans-load--------------------------------------------------------------
# Load autos dataset
data(autos)

# Display structure
cat("Dimensions:", dim(autos), "\n")
cat("Variables:", paste(head(names(autos), 10), collapse = ", "), "...\n")

# Summary
summary(autos[, 1:6])

## ----kmeans-fit---------------------------------------------------------------
# Initialize KmeansVariables
km <- KmeansVariables$new(
  n_clusters = 5,
  n_init = 10,
  distance_metric = "r_squared",
  random_state = 42
)

# Fit
km$fit(autos)

# Print
print(km)

## ----kmeans-results-----------------------------------------------------------
# Cluster assignments
head(km$labels, 10)

# Cluster sizes
km$cluster_sizes

# Total inertia
cat("Total inertia:", round(km$inertia, 4), "\n")

## ----kmeans-elbow-------------------------------------------------------------
km$plot_elbow()

## ----kmeans-projection--------------------------------------------------------
km$plot_projection(show_labels = TRUE)

## ----kmeans-contributions-----------------------------------------------------
km$plot_contributions(top_n = 5)

## ----kmeans-summary-----------------------------------------------------------
results_km <- km$summary(print_summary = TRUE)

## ----kmeans-clust-summary-----------------------------------------------------
results_km$clust_summary

## ----kmeans-members-----------------------------------------------------------
head(results_km$clust_members, 10)

## ----kmeans-centroids---------------------------------------------------------
round(results_km$centroids_correlations, 3)

## ----modclust-load------------------------------------------------------------
# Load vote dataset
data(vote)

# Display
cat("Dimensions:", dim(vote), "\n")
cat("Variables:", paste(names(vote), collapse = ", "), "\n")

# Summary
summary(vote)

## ----modclust-split-----------------------------------------------------------
# Split for prediction example
set.seed(42)
train_idx <- sample(nrow(vote), size = round(0.8 * nrow(vote)))
vote_train <- vote[train_idx, ]
vote_test <- vote[-train_idx, ]

cat("Training:", nrow(vote_train), "observations\n")
cat("Test:", nrow(vote_test), "observations\n")

## ----modclust-fit-------------------------------------------------------------
# Initialize ModCluster
mc <- ModCluster$new(
  method = "hclust",
  n_dimensions = 2,
  hclust_method = "average"
)

# Fit
mc$fit(vote_train)

# Print
print(mc)

## ----modclust-heights---------------------------------------------------------
mc$plot_heights()

## ----modclust-cut-------------------------------------------------------------
mc$cut_tree(3)
print(mc)

## ----modclust-dendro----------------------------------------------------------
mc$plot_dendrogram()

## ----modclust-mca-------------------------------------------------------------
mc$plot_mca()

## ----modclust-summary---------------------------------------------------------
summary(mc)

## ----session-info-------------------------------------------------------------
sessionInfo()

