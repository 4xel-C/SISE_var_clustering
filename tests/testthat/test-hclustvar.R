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


# ==============================================================================
# Test: Cramer's calculation method
# ==============================================================================


test_that("private$cramer_v() computes expected Cramer's V values", {
  # Create an instance
  obj <- HClustVar$new(vartype = "qual")

  # Access the private method
  cramer_v <- obj$.__enclos_env__$private$cramer_v

  # Simple association: identical variables → V = 1
  x <- c("A", "B", "A", "B")
  y <- c("A", "B", "A", "B")
  expect_equal(round(cramer_v(x, y), 5), 1)

  # Independent variables → low V
  x <- c("A", "A", "B", "B")
  y <- c("X", "Y", "X", "Y")
  v <- cramer_v(x, y)
  expect_true(v >= 0 && v <= 1)

  # Different number of levels
  x <- factor(c("A", "A", "B", "C"))
  y <- factor(c("X", "Y", "X", "Y"))
  v2 <- cramer_v(x, y)
  expect_true(is.numeric(v2))
  expect_true(v2 >= 0 && v2 <= 1)
})

# ==============================================================================
# Test: Cramer's matrix generation method
# ==============================================================================

test_that("private$cramer_matrix() returns valid symmetric matrix", {
  # Create an instance
  obj <- HClustVar$new(vartype = "qual")

  # Access the private method
  cramer_matrix <- obj$.__enclos_env__$private$cramer_matrix

  # Create dummy data
  df <- data.frame(
    A = c("A", "A", "B", "B"),
    B = c("X", "Y", "X", "Y"),
    C = c("M", "M", "N", "N")
  )

  # Compute matrix
  mat <- cramer_matrix(df)

  # Expect a square symmetric matrix
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), ncol(mat))
  expect_true(all(rownames(mat) == colnames(mat)))
  expect_equal(mat, t(mat))  # Symmetry check

  # Diagonal values should be 1
  expect_true(all(diag(mat) == rep(1, ncol(df))))

  # Off-diagonal values should be between 0 and 1
  off_diag <- mat[upper.tri(mat)]
  expect_true(all(off_diag >= 0 & off_diag <= 1))
})

test_that("cramer_matrix() integrates cramer_v correctly", {
  obj <- HClustVar$new(vartype = "qual")

  cramer_v <- obj$.__enclos_env__$private$cramer_v
  cramer_matrix <- obj$.__enclos_env__$private$cramer_matrix

  df <- data.frame(
    A = c("A", "A", "B", "B"),
    B = c("A", "B", "A", "B")
  )

  mat <- cramer_matrix(df)

  # Expected diagonal = 1
  expect_equal(diag(mat), c(A = 1, B = 1))

  # Expected Cramer's V off-diagonal equals direct computation
  v_direct <- cramer_v(df$A, df$B)
  expect_equal(mat["A", "B"], v_direct)
  expect_equal(mat["B", "A"], v_direct)
})

# ==============================================================================
# TEST: Quantile discretisation method
# ==============================================================================

# Creating minimal instance of the class.
hclust <- HClustVar$new(vartype = "mixed")

# Exemple dataframe.
df <- data.frame(
  var1 = 1:10,
  var2 = c(10, 9, 1, 2, 3, 8, 6, 7, 4, 5),
  var3 = c(rep(1, 5), rep(2, 5))
)

# Access to private method for testing.
quantile_discretisation <- hclust$.__enclos_env__$private$quantile_discretisation

test_that("quantile_discretisation transform quantitative into qualitative", {
  df_disc <- quantile_discretisation(df, quanti_index = c(1, 2), n_groups = 4)

  # Quantitatives should become factors.
  expect_true(is.factor(df_disc$var1))
  expect_true(is.factor(df_disc$var2))

  # We should have Q1 to Q4 labels.
  expect_equal(levels(df_disc$var1), c("Q1", "Q2", "Q3", "Q4"))

  # Last column unchanged.
  expect_false(is.factor(df_disc$var3))
  expect_equal(df_disc$var3, df$var3)
})

test_that("quantile_discretisation n_groups behave correctly", {
  df_disc <- quantile_discretisation(df, quanti_index = c("var1", "var2"), n_groups = 3)
  expect_equal(levels(df_disc$var1), c("Q1", "Q2", "Q3"))

  df_disc <- quantile_discretisation(df, quanti_index = c("var1", "var2"), n_groups = 2)
  expect_equal(levels(df_disc$var1), c("Q1", "Q2"))

  df_disc <- quantile_discretisation(df, quanti_index = c("var1", "var2"), n_groups = 1)
  expect_equal(levels(df_disc$var1), c("Q1"))
})

test_that("quantile_discretisation keep dataframe size", {
  df_disc <- quantile_discretisation(df, quanti_index = c("var1", "var2"), n_groups = 4)
  expect_equal(dim(df_disc), dim(df))
})

test_that("quantile_discretisation manage missing values.", {
  df_na <- df
  df_na$var1[c(3, 7)] <- NA
  expect_no_error({
    df_disc <- quantile_discretisation(df_na, quanti_index = "var1", n_groups = 4)
  })
})

test_that("quantile_discretisation manage only specified group.", {
  df_disc <- quantile_discretisation(df, quanti_index = "var1", n_groups = 1)
  expect_equal(levels(df_disc$var1), "Q1")
})

test_that("quantile_discretisation return identical dataframe if no column.", {
  df_disc <- quantile_discretisation(df, quanti_index = integer(0), n_groups = 4)
  expect_equal(df_disc, df)
})


test_that("quantile_discretisation return the correct labeled values.", {
  df_disc <- quantile_discretisation(df, quanti_index = c(1, 2), n_groups = 4)

  expect_equal(df_disc, data.frame(
    var1 = factor(c(rep("Q1", 3), rep("Q2", 2), rep("Q3", 2), rep("Q4", 3))),
    var2 = factor(c("Q4", "Q4", "Q1", "Q1", "Q1", "Q4", "Q3", "Q3", "Q2", "Q2")),
    var3 = c(rep(1, 5), rep(2, 5))
  ))
})
