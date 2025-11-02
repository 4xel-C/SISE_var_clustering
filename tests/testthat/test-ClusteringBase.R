# ==============================================================================
# Unit Tests for ClusteringBase (Abstract Class)
# ==============================================================================

# ==============================================================================
# Helper: Create a concrete test class that inherits from ClusteringBase
# ==============================================================================

# Dummy child class
DummyClustering <- R6::R6Class(
  "DummyClustering",
  inherit = ClusteringBase,
  public = list(
    fit = function() { TRUE },
    predict = function() { TRUE },
    summary = function() { TRUE },
    print = function() { TRUE }
  )
)

# Test datasets
df_mixed <- data.frame(
  num1 = 1:5,
  num2 = c(5,4,3,2,1),
  cat1 = factor(c("a","b","a","b","a"))
)

df_quanti <- data.frame(num1 = 1:5, num2 = 6:10)
df_quali  <- data.frame(cat1 = factor(c("a","b","a","b","a")), cat2 = factor(c("x","y","x","y","x")))


# ==============================================================================
# Test Suite: Abstract Class Prevention
# ==============================================================================

test_that("ClusteringBase cannot be instantiated directly", {
  expect_error(
    ClusteringBase$new(),
    "ClusteringBase is an abstract class and shouldn't be directly instanciated."
  )
})

test_that("Child class can be instantiated", {
  obj <- DummyClustering$new()
  expect_s3_class(obj, "DummyClustering")
})


# ==============================================================================
# Test: Data Structure Validation
# ==============================================================================

test_that("load_and_check_data works and detects types", {
  obj <- DummyClustering$new()
  obj$load_and_check_data(df_mixed)

  expect_equal(obj$quanti_indices, c(1,2))
  expect_equal(obj$quali_indices, 3)
  expect_equal(nrow(obj$get_quanti_data()), nrow(df_mixed))
  expect_equal(ncol(obj$get_quali_data()), 1)
})


test_that("n_clusters setter validates correctly", {
  obj <- DummyClustering$new()
  obj$load_and_check_data(df_quanti)

  obj$n_clusters <- 2
  expect_equal(obj$n_clusters, 2)
  expect_error(obj$n_clusters <- 0)
  expect_error(obj$n_clusters <- 10) # > number of variables
  expect_error(obj$n_clusters <- "a")
})


test_that("labels setter validates length", {
  obj <- DummyClustering$new()
  obj$load_and_check_data(df_quanti)

  obj$labels <- 1:5
  expect_equal(obj$labels, 1:5)
  expect_error(obj$labels <- 1:3)
})


test_that("fitted setter works", {
  obj <- DummyClustering$new()
  obj$fitted <- TRUE
  expect_true(obj$fitted)
  expect_error(obj$fitted <- 1)
})


test_that("check_missing_values stops if NA present", {
  obj <- DummyClustering$new()
  df_na <- df_quanti
  df_na[1,1] <- NA
  expect_error(obj$load_and_check_data(df_na))
})

test_that("1 column dataframe stops", {
  obj <- DummyClustering$new()
  df_test <- df_quanti[1]
  expect_error(obj$load_and_check_data(df_test), "At least 2 variables are required for clustering, got: 1")
})

test_that("1 row dataframe stops", {
  obj <- DummyClustering$new()
  df_test <- df_quanti[1,]
  expect_error(obj$load_and_check_data(df_test), "At least 2 observations are required, got: 1")
})

test_that("Dataframe with no name got default column names.", {
  obj <- DummyClustering$new()
  df_test <- unname(df_quanti)
  expect_warning(obj$load_and_check_data(df_test), "No column names provided. Using default names: V1, V2, ...")
  expect_equal(colnames(obj$data), c("V1", "V2"))
})
