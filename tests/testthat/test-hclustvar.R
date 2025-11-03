# ==============================================================================
# Unit Tests for HClustVar
# ==============================================================================

# ==============================================================================
# Constructor.
# ==============================================================================
test_that("Constructor sets default values correctly", {
  obj <- HClustVar$new()
  expect_equal(obj$vartype, "mixed")
  expect_null(obj$dist.metric)
})

test_that("Constructor sets quant type with default dist.metric", {
  obj <- HClustVar$new(vartype = "quant")
  expect_equal(obj$vartype, "quant")
  expect_equal(obj$dist.metric, "rsquare")
})

test_that("Constructor sets quant type with custom dist.metric", {
  obj <- HClustVar$new(vartype = "quant", dist.metric = "r")
  expect_equal(obj$dist.metric, "r")
})

test_that("Constructor warns when dist.metric is ignored for qual/mixed", {
  expect_warning(obj1 <- HClustVar$new(vartype = "qual", dist.metric = "r"))
  expect_equal(obj1$vartype, "qual")
  expect_equal(obj1$dist.metric, "r")

  expect_warning(obj2 <- HClustVar$new(vartype = "mixed", dist.metric = "rsquare"))
  expect_equal(obj2$vartype, "mixed")
})

test_that("check_input() errors for invalid vartype or dist.metric", {
  obj <- HClustVar$new()

  # Vérifie vartype invalide
  expect_error(
    obj$initialize(vartype = "wrong"),
    "Parameter 'vartype' has invalid value"
  )

  # Vérifie dist.metric invalide pour quant
  expect_error(
    obj$initialize(vartype = "quant", dist.metric = "invalid"),
    "Parameter 'dist.metric' has invalid value"
  )
})
