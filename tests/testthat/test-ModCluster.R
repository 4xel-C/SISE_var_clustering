# Test unitaires pour ModCluster
# Fichier: testthat/test-ModCluster.R

library(testthat)
library(R6)
library(FactoMineR)

# ==============================================================================
# DONNÉES DE TEST
# ==============================================================================

# Fonction pour créer des données de test
create_test_data <- function(n = 100, seed = 123) {
  set.seed(seed)
  data.frame(
    color = factor(sample(c("red", "blue", "green"), n, replace = TRUE)),
    size = factor(sample(c("small", "medium", "large"), n, replace = TRUE)),
    shape = factor(sample(c("circle", "square", "triangle"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
}

# Données pour les tests
test_data <- create_test_data(100)
test_data_small <- create_test_data(30, seed = 456)

# ==============================================================================
# TESTS D'INITIALISATION
# ==============================================================================

test_that("ModCluster initialisation avec paramètres par défaut", {
  model <- ModCluster$new()

  expect_s3_class(model, "ModCluster")
  expect_s3_class(model, "R6")
  expect_false(model$fitted)
})


test_that("ModCluster initialisation avec différentes méthodes hclust", {
  methods <- c("average", "single", "complete")

  for (method in methods) {
    model <- ModCluster$new(hclust_method = method)
    expect_s3_class(model, "ModCluster")
  }
})

test_that("ModCluster rejette les méthodes hclust invalides", {
  expect_error(
    ModCluster$new(hclust_method = "invalid_method"),
    "Invalid hclust_method"
  )
})

test_that("ModCluster rejette les méthodes de clustering non supportées", {
  expect_error(
    ModCluster$new(method = "kmeans"),
    "Only 'hclust' is supported"
  )
})

# ==============================================================================
# TESTS DE FIT
# ==============================================================================

test_that("fit() fonctionne avec des données valides", {
  model <- ModCluster$new()

  expect_message(model$fit(test_data), "completed successfully")
  expect_true(model$fitted)
})

test_that("fit() crée la matrice disjonctive correctement", {
  model <- ModCluster$new()
  model$fit(test_data)

  expect_false(is.null(model$modality_names))
  expect_equal(length(model$modality_names), 9) # 3 + 3 + 3 niveaux

  # Vérifier le format des noms de modalités
  expect_true(all(grepl("\\.", model$modality_names)))
})

test_that("fit() calcule les fréquences des modalités", {
  model <- ModCluster$new()
  model$fit(test_data)

  freqs <- model$modality_frequencies
  expect_false(is.null(freqs))
  expect_equal(length(freqs), 9)
  expect_true(all(freqs > 0))
  expect_equal(sum(freqs), nrow(test_data) * 3) # 3 variables
})

test_that("fit() calcule la matrice de distances Dice", {
  model <- ModCluster$new()
  model$fit(test_data)

  dice_mat <- model$get_dice_matrix()
  expect_false(is.null(dice_mat))
  expect_equal(nrow(dice_mat), 9)
  expect_equal(ncol(dice_mat), 9)
  expect_true(isSymmetric(dice_mat))
  # FIX 1: Comparer les valeurs numériques sans les noms
  expect_equal(as.numeric(diag(dice_mat)), rep(0, 9))
})

test_that("fit() calcule l'ACM pour la visualisation", {
  model <- ModCluster$new(, n_dimensions = 5)
  model$fit(test_data)

  mca <- model$mca_result
  expect_false(is.null(mca))
  expect_s3_class(mca, "MCA")
})

test_that("fit() échoue avec des données non catégorielles", {
  numeric_data <- data.frame(x = rnorm(100), y = rnorm(100))
  model <- ModCluster$new()

  # FIX 2: Utiliser le message d'erreur réel
  expect_error(
    model$fit(numeric_data),
    "Not enough qualitative columns"
  )
})

# ==============================================================================
# TESTS DE CUT_TREE
# ==============================================================================

test_that("cut_tree() fonctionne après fit()", {
  model <- ModCluster$new()
  model$fit(test_data)

  expect_message(model$cut_tree(3), "Tree cut into 3 clusters")
  expect_equal(length(unique(model$modality_labels)), 3)
})

test_that("cut_tree() modifie le nombre de clusters", {
  model <- ModCluster$new()
  model$fit(test_data)

  model$cut_tree(2)
  expect_equal(length(unique(model$modality_labels)), 2)

  model$cut_tree(4)
  expect_equal(length(unique(model$modality_labels)), 4)
})

test_that("cut_tree() échoue si le modèle n'est pas fitted", {
  model <- ModCluster$new()

  expect_error(
    model$cut_tree(3),
    "Model must be fitted"
  )
})

# ==============================================================================
# TESTS DE PREDICT
# ==============================================================================

test_that("predict() fonctionne avec de nouvelles données", {
  model <- ModCluster$new()
  model$fit(test_data)
  model$cut_tree(3)

  new_data <- create_test_data(100, seed = 999)
  predictions <- model$predict(new_data)$prediction

  expect_true(all(predictions %in% 1:3))
})

test_that("predict() échoue si le modèle n'est pas fitted", {
  model <- ModCluster$new()
  new_data <- create_test_data(100)

  expect_error(
    model$predict(new_data),
    "Model must be fitted"
  )
})

test_that("predict() échoue si les clusters ne sont pas définis", {
  model <- ModCluster$new()
  model$fit(test_data)
  new_data <- create_test_data(100)

  expect_error(
    model$predict(new_data),
    "Tree must be cutted first"
  )
})


# ==============================================================================
# TESTS DE PROJECT_NEW_MODALITIES
# ==============================================================================

test_that("project_new_modalities() projette correctement", {
  model <- ModCluster$new()
  model$fit(test_data)

  new_data <- create_test_data(5, seed = 888)
  projected <- model$project_new_modalities(new_data)

  expect_true(is.matrix(projected))
  expect_equal(nrow(projected), 5)
  expect_true(all(grepl("^Dim\\.", colnames(projected))))
})

test_that("project_new_modalities() respecte les axes spécifiés", {
  model <- ModCluster$new(, n_dimensions = 5)
  model$fit(test_data)

  new_data <- create_test_data(5, seed = 888)
  projected <- model$project_new_modalities(new_data, axes = c(1, 3))

  expect_equal(ncol(projected), 2)
  expect_equal(colnames(projected), c("Dim.1", "Dim.3"))
})

test_that("project_new_modalities() échoue si non fitted", {
  model <- ModCluster$new()
  new_data <- create_test_data(5)

  expect_error(
    model$project_new_modalities(new_data),
    "Model must be fitted"
  )
})

# ==============================================================================
# TESTS DES MÉTHODES D'INFORMATION
# ==============================================================================

test_that("summary() retourne les bonnes informations", {
  model <- ModCluster$new()
  model$fit(test_data)
  model$cut_tree(3)

  summ <- model$summary()

  expect_type(summ, "list")
  expect_true("clust_summary" %in% names(summ))
  expect_true("clust_members" %in% names(summ))
  expect_equal(nrow(summ$clust_summary), 3)
  expect_equal(nrow(summ$clust_members), 9)
})

test_that("get_cluster_table() retourne un tableau correct", {
  model <- ModCluster$new()
  model$fit(test_data)

  model$cut_tree(4)

  table <- model$get_cluster_table()

  expect_s3_class(table, "data.frame")
  expect_equal(nrow(table), 9)
  expect_true(all(c("Modality", "Cluster", "Frequency") %in% names(table)))
})

test_that("get_dice_matrix() retourne la matrice de distances", {
  model <- ModCluster$new()
  model$fit(test_data)

  dice_mat <- model$get_dice_matrix()

  expect_true(is.matrix(dice_mat))
  expect_equal(dim(dice_mat), c(9, 9))
  expect_true(isSymmetric(dice_mat))
})

# ==============================================================================
# TESTS DES ACTIVE BINDINGS
# ==============================================================================

test_that("Active bindings retournent les bonnes valeurs", {
  model <- ModCluster$new()
  model$fit(test_data)

  expect_true(is.null(model$modality_labels))
  expect_false(is.null(model$modality_names))
  expect_false(is.null(model$modality_frequencies))
  expect_false(is.null(model$mca_result))
})

# ==============================================================================
# TESTS D'INTÉGRATION
# ==============================================================================

test_that("Workflow avec cut_tree() après fit", {
  model <- ModCluster$new()
  model$fit(test_data)

  # Coupe l'arbre
  model$cut_tree(4)

  # FIX 4: Vérifier que le modèle a bien 4 clusters
  expect_equal(length(unique(model$modality_labels)), 4)

  # Maintenant ça fonctionne
  pred <- model$predict(test_data)

  # pas toutes les observations
  # utilisent nécessairement tous les clusters
  n_predicted_clusters <- length(unique(pred$prediction))
  expect_true(n_predicted_clusters >= 1 && n_predicted_clusters <= 4,
              info = paste("Expected 1-4 clusters, got", n_predicted_clusters))
})

# ==============================================================================
# TESTS DE ROBUSTESSE
# ==============================================================================

test_that("Gère les données avec peu d'observations", {
  small_data <- create_test_data(10, seed = 222)
  model <- ModCluster$new()

  expect_message(model$fit(small_data))
  expect_true(model$fitted)
})

test_that("Gère différents nombres de clusters", {
  model <- ModCluster$new()
  model$fit(test_data)

  for (k in 2:6) {
    model$cut_tree(k)
    expect_equal(length(unique(model$modality_labels)), k)
  }
})

cat("\n✓ Tous les tests sont définis et prêts à être exécutés!\n")
