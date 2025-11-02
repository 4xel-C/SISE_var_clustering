# ==============================================================================
# TESTS UNITAIRES POUR ClusteringBase
# ==============================================================================
# File: tests/testthat/test-ClusteringBase.R

library(testthat)
library(R6)

# ------------------------------------------------------------------------------
# HELPER: Create a concrete implementation for testing abstract class
# ------------------------------------------------------------------------------

DummyClustering <- R6::R6Class(
  "DummyClustering",
  inherit = ClusteringBase,
  public = list(
    fit = function() { private$fitted <- TRUE },
    predict = function() { return(rep(1, nrow(private$.data))) },
    summary = function() { return("Dummy summary") },
    print = function() { cat("<DummyClustering object>\n") }
  )
)

# ------------------------------------------------------------------------------
# TEST: ABSTRACT CLASS INSTANTIATION
# ------------------------------------------------------------------------------

test_that("ClusteringBase cannot be instantiated directly", {
  expect_error(
    ClusteringBase$new(),
    "ClusteringBase is an abstract class"
  )
})

test_that("Concrete child class can be instantiated", {
  data_test <- data.frame(x = 1:10, y = 11:20)

  expect_silent({
    model <- DummyClustering$new(data = data_test, n_clusters = 2)
  })

  expect_s3_class(model, "DummyClustering")
  expect_s3_class(model, "ClusteringBase")
})

# ------------------------------------------------------------------------------
# TEST 2: DATA VALIDATION
# ------------------------------------------------------------------------------

test_that("Data must be data.frame or matrix", {
  expect_error(
    DummyClustering$new(data = list(a = 1:5, b = 6:10), n_clusters = 2),
    "'data' must be a data.frame or matrix"
  )

  expect_error(
    DummyClustering$new(data = c(1, 2, 3, 4, 5), n_clusters = 2),
    "'data' must be a data.frame or matrix"
  )

  expect_error(
    DummyClustering$new(data = "not_a_dataframe", n_clusters = 2),
    "'data' must be a data.frame or matrix"
  )
})

test_that("Matrix is converted to data.frame", {
  mat <- matrix(1:20, nrow = 10, ncol = 2)
  colnames(mat) <- c("V1", "V2")

  model <- DummyClustering$new(data = mat, n_clusters = 2)

  expect_true(is.data.frame(model$data))
  expect_equal(ncol(model$data), 2)
  expect_equal(nrow(model$data), 10)
})


test_that("Data must have at least 2 variables", {
  expect_error(
    DummyClustering$new(data = data.frame(x = 1:10), n_clusters = 1),
    "At least 2 variables are required for clustering, got: 1"
  )
})



test_that("Data must have at least 2 observations", {
  expect_error(
    DummyClustering$new(data = data.frame(x = 1, y = 2), n_clusters = 2),
    "At least 2 observations are required"
  )
})


test_that("Valid data.frame is accepted", {
  data_test <- data.frame(
    var1 = rnorm(20),
    var2 = rnorm(20),
    var3 = rnorm(20)
  )

  expect_silent({
    model <- DummyClustering$new(data = data_test, n_clusters = 2)
  })

  expect_equal(model$data, data_test)
})




# ------------------------------------------------------------------------------
# TEST 3: MISSING VALUES
# ------------------------------------------------------------------------------

test_that("Missing values trigger an error", {
  data_with_na <- data.frame(
    x = c(1, 2, NA, 4, 5),
    y = c(10, 20, 30, NA, 50)
  )

  expect_error(
    DummyClustering$new(data = data_with_na, n_clusters = 2),
    "Your data contains NA values"
  )
})

test_that("Data without NA is accepted", {
  data_no_na <- data.frame(
    x = 1:10,
    y = 11:20,
    z = 21:30
  )

  expect_silent({
    model <- DummyClustering$new(data = data_no_na, n_clusters = 2)
  })
})


# ------------------------------------------------------------------------------
# TEST 4: VARIABLE TYPE DETECTION
# ------------------------------------------------------------------------------


test_that("Qualitative variables are detected correctly", {
  data_quali <- data.frame(
    var1 = factor(sample(c("A", "B", "C"), 20, replace = TRUE)),
    var2 = factor(sample(c("X", "Y"), 20, replace = TRUE)),
    var3 = sample(c("Low", "High"), 20, replace = TRUE)  # character
  )

  model <- DummyClustering$new(data = data_quali, n_clusters = 2)

  expect_equal(length(model$quanti_indices), 0)
  expect_equal(length(model$quali_indices), 3)
  expect_equal(model$quali_indices, 1:3)
})

test_that("Quantitative variables are detected correctly", {
  data_quali <- data.frame(
    var1 = sample(c(1, 2, 3), 20, replace = TRUE),
    var2 = sample(c(6, 8), 20, replace = TRUE),
    var3 = sample(c(3, 4), 20, replace = TRUE)
  )

  model <- DummyClustering$new(data = data_quali, n_clusters = 2)

  expect_equal(length(model$quanti_indices), 3)
  expect_equal(length(model$quali_indices), 0)
  expect_equal(model$quali_indices, integer(0))
  expect_equal(model$quanti_indices, 1:3)
})

test_that("Variable type detection works", {
  d <- data.frame(x = 1:5, y = factor(c("a", "b", "a", "b", "a")), z = 6:10, a = c("one", "two", "three", "four", "five"))
  obj <- DummyClustering$new(d)
  expect_equal(obj$quanti_indices, c(1, 3))
  expect_equal(obj$quali_indices, c(2, 4))
})


# ------------------------------------------------------------------------------
# TEST 4: METHODS
# ------------------------------------------------------------------------------

test_that("get_quanti_data and get_quali_data behave correctly", {
  d <- data.frame(x = 1:5, y = factor(c("a", "b", "a", "b", "a")), z = 6:10, a = c("one", "two", "three", "four", "five"))
  obj <- DummyClustering$new(d)
  expect_true(is.data.frame(obj$get_quanti_data()))
  expect_true(is.data.frame(obj$get_quali_data()))
  expect_equal(ncol(obj$get_quanti_data()), 2)
  expect_equal(ncol(obj$get_quali_data()), 2)
  expect_equal(obj$get_quanti_data(), data.frame(x = 1:5, z = 6:10))
  expect_equal(obj$get_quali_data(), data.frame(y = factor(c("a", "b", "a", "b", "a")), a = c("one", "two", "three", "four", "five")))
})


test_that("get_quanti_data and get_quali_data always return a dataframe (even when length == 1).", {
  d <- data.frame(x = 1:5, a = c("one", "two", "three", "four", "five"))
  obj <- DummyClustering$new(d)
  expect_true(is.data.frame(obj$get_quanti_data()))
  expect_true(is.data.frame(obj$get_quali_data()))
  expect_equal(ncol(obj$get_quanti_data()), 1)
  expect_equal(ncol(obj$get_quali_data()), 1)
})


# ------------------------------------------------------------------------------
# TEST 5: Getter/Setters
# ------------------------------------------------------------------------------

test_that("data getter work correctly", {
  df <- data.frame(x = 1:3, y = 4:6)
  cb <- DummyClustering$new(df)
  expect_equal(cb$data, df)
})

test_that("n_clusters setter validates input", {
  cb <- DummyClustering$new()
  cb$n_clusters <- 3
  expect_equal(cb$n_clusters, 3)

  expect_error(cb$n_clusters <- -1, "'n_clusters' must be a positive integer scalar")
  expect_error(cb$n_clusters <- "a", "'n_clusters' must be a positive integer scalar")
})

test_that("labels setter validates length", {
  df <- data.frame(x = 1:3, y = 4:6)
  cb <- DummyClustering$new(df)

  cb$labels <- c(1, 2, 3)
  expect_equal(cb$labels, c(1, 2, 3))

  expect_error(cb$labels <- c(1, 2), "'labels' length must match number of observations.")
})

test_that("fitted setter only accepts logical of length 1", {
  cb <- DummyClustering$new()
  cb$fitted <- TRUE
  expect_true(cb$fitted)

  expect_error(cb$fitted <- c(TRUE, FALSE), "single logical value")
  expect_error(cb$fitted <- "TRUE", "logical value")
})

test_that("quanti_indices and quali_indices getters return expected defaults", {
  cb <- DummyClustering$new()
  expect_null(cb$quanti_indices)
  expect_null(cb$quali_indices)
})
