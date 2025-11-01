# ==============================================================================
# TESTS UNITAIRES POUR ClusteringBase
# ==============================================================================
# File: tests/testthat/test-ClusteringBase.R

library(testthat)
library(R6)

# ------------------------------------------------------------------------------
# HELPER: Create a concrete implementation for testing abstract class
# ------------------------------------------------------------------------------

# Concrete implementation of ClusteringBase for testing purposes
TestClustering <- R6::R6Class(
  "TestClustering",
  inherit = ClusteringBase,

  public = list(

    fit = function() {
      private$.fitted <- TRUE
      private$.labels <- sample(1:private$.n_clusters,
                                ncol(private$.data),
                                replace = TRUE)
      names(private$.labels) <- colnames(private$.data)
      invisible(self)
    },

    predict = function(newdata = NULL) {
      if (is.null(newdata)) {
        return(private$.labels)
      }
      return(sample(1:private$.n_clusters, ncol(newdata), replace = TRUE))
    },

    summary = function() {
      cat("TestClustering Summary\n")
      invisible(self)
    },

    print = function() {
      cat("<TestClustering>\n")
      invisible(self)
    }
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
    model <- TestClustering$new(data = data_test, n_clusters = 2)
  })

  expect_s3_class(model, "TestClustering")
  expect_s3_class(model, "ClusteringBase")
})

# ------------------------------------------------------------------------------
# TEST SUITE 2: DATA VALIDATION
# ------------------------------------------------------------------------------

test_that("Data must be data.frame or matrix", {
  expect_error(
    TestClustering$new(data = list(a = 1:5, b = 6:10), n_clusters = 2),
    "'data' must be a data.frame or matrix"
  )

  expect_error(
    TestClustering$new(data = c(1, 2, 3, 4, 5), n_clusters = 2),
    "'data' must be a data.frame or matrix"
  )

  expect_error(
    TestClustering$new(data = "not_a_dataframe", n_clusters = 2),
    "'data' must be a data.frame or matrix"
  )
})

test_that("Matrix is converted to data.frame", {
  mat <- matrix(1:20, nrow = 10, ncol = 2)
  colnames(mat) <- c("V1", "V2")

  model <- TestClustering$new(data = mat, n_clusters = 2)

  expect_true(is.data.frame(model$data))
  expect_equal(ncol(model$data), 2)
  expect_equal(nrow(model$data), 10)
})


test_that("Data must have at least 2 variables", {
  expect_error(
    TestClustering$new(data = data.frame(x = 1:10), n_clusters = 1),
    "At least 2 variables are required for clustering, got: 1"
  )
})

test_that("Data must have at least 2 observations", {
  expect_error(
    TestClustering$new(data = data.frame(x = 1, y = 2), n_clusters = 2),
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
    model <- TestClustering$new(data = data_test, n_clusters = 2)
  })

  expect_equal(model$data, data_test)
})

# ------------------------------------------------------------------------------
# TEST SUITE 4: MISSING VALUES
# ------------------------------------------------------------------------------

test_that("Missing values trigger an error", {
  data_with_na <- data.frame(
    x = c(1, 2, NA, 4, 5),
    y = c(10, 20, 30, NA, 50)
  )

  expect_error(
    TestClustering$new(data = data_with_na, n_clusters = 2),
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
    model <- TestClustering$new(data = data_no_na, n_clusters = 2)
  })
})


# ------------------------------------------------------------------------------
# TEST SUITE 5: VARIABLE TYPE DETECTION
# ------------------------------------------------------------------------------


test_that("Qualitative variables are detected correctly", {
  data_quali <- data.frame(
    var1 = factor(sample(c("A", "B", "C"), 20, replace = TRUE)),
    var2 = factor(sample(c("X", "Y"), 20, replace = TRUE)),
    var3 = sample(c("Low", "High"), 20, replace = TRUE)  # character
  )

  model <- TestClustering$new(data = data_quali, n_clusters = 2)

  expect_equal(length(model$quanti_indices), 0)
  expect_equal(length(model$quali_indices), 3)
  expect_equal(model$quali_indices, 1:3)
})

