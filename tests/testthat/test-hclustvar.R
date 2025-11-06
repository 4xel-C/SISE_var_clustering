# ==============================================================================
# Unit Tests for HClustVar
# ==============================================================================

# ==============================================================================
# Constructor and fit method.
# ==============================================================================
test_that("Constructor sets default values correctly", {
  obj <- HClustVar$new()
  expect_equal(obj$vartype, "auto")
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


test_that("HClustVar auto detection with quantitative data", {
  df_quanti <- data.frame(
    var1 = rnorm(10),
    var2 = rnorm(10),
    var3 = rnorm(10)
  )

  hc_auto_quanti <- HClustVar$new(vartype = "auto")
  expect_silent(hc_auto_quanti$fit(df_quanti))

  expect_equal(hc_auto_quanti$vartype, "quant")
  expect_equal(hc_auto_quanti$dist.metric, "rsquare")
  expect_equal(hc_auto_quanti$fitted, TRUE)
  expect_s3_class(hc_auto_quanti$dist.matrix, "dist")
})


test_that("HClustVar auto detection with qualitative data", {
  df_quali <- data.frame(
    color = factor(sample(c("red", "blue", "green"), 10, TRUE)),
    shape = factor(sample(c("circle", "square"), 10, TRUE))
  )

  hc_auto_quali <- HClustVar$new(vartype = "auto")
  expect_silent(hc_auto_quali$fit(df_quali))

  expect_equal(hc_auto_quali$vartype, "qual")
  expect_null(hc_auto_quali$dist.metric)
  expect_equal(hc_auto_quali$fitted, TRUE)
  expect_s3_class(hc_auto_quali$dist.matrix, "dist")
})


test_that("HClustVar auto detection with mixed data", {
  df_mixed <- data.frame(
    age = rnorm(10),
    income = rnorm(10),
    gender = factor(sample(c("M", "F"), 10, TRUE))
  )

  hc_auto_mixed <- HClustVar$new(vartype = "auto")
  expect_silent(hc_auto_mixed$fit(df_mixed))

  expect_equal(hc_auto_mixed$vartype, "mixed")
  expect_null(hc_auto_mixed$dist.metric)
  expect_s3_class(hc_auto_mixed$dist.matrix, "dist")
})


test_that("HClustVar handles invalid metric gracefully", {
  df_quali <- data.frame(
    color = factor(sample(c("red", "blue", "green"), 10, TRUE)),
    shape = factor(sample(c("circle", "square"), 10, TRUE))
  )

  hc_invalid <- HClustVar$new(vartype = "auto", dist.metric = "r")

  expect_warning(
    hc_invalid$fit(df_quali),
    regexp = "Dataframe has only qualitative values, param. dist.metric is ignored."
  )
  expect_equal(hc_invalid$fitted, TRUE)
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

# Dumified dataframe
df01 <- data.frame(
  var1 = 1:10,
  var2 = c(10, 9, 1, 2, 3, 8, 6, 7, 4, 5),
  var3 = c(0, 1, 0, 0, 1, 0, 1, 1, 1, 0)
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

test_that("quantile_discretisation change dummified columns to factor without changing values.", {
  df_disc <- quantile_discretisation(df01, quanti_index = c("var1", "var2", "var3"), n_groups = 4)
  expect_equal(df_disc$var3, as.factor(df01$var3))
})

# ==============================================================================
# Test cut-tree method.
# ==============================================================================

test_that("cut_tree returns correct cluster labels with k", {
  # Create and HClustVar object for tests.
  obj <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
  data_quant <- data.frame(
    var1 = rnorm(10),
    var2 = rnorm(10),
    var3 = rnorm(10)
  )
  obj$fit(data_quant)

  k <- 2
  labels <- obj$cut_tree(k = k)

  # Check length to be equal to the number of columns or the correct subset.
  expect_equal(length(labels), ncol(data_quant))

  # Check that all labels are included in the correct range.
  expect_true(all(labels %in% 1:k))
})

test_that("cut_tree returns correct cluster labels with h", {
  obj <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
  data_quant <- data.frame(
    var1 = rnorm(10),
    var2 = rnorm(10),
    var3 = rnorm(10)
  )
  obj$fit(data_quant)

  # Cut at 0.5 height.
  h <- 0.5
  labels <- obj$cut_tree(h = h)

  # Check length
  expect_equal(length(labels), ncol(data_quant))

  # Check all labels are positive.
  expect_true(all(labels > 0))
})

test_that("cut_tree prioritizes k over h if both provided", {
  obj <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
  data_quant <- data.frame(
    var1 = rnorm(10),
    var2 = rnorm(10),
    var3 = rnorm(10)
  )
  obj$fit(data_quant)

  k <- 2
  h <- 0.1
  labels <- obj$cut_tree(k = k, h = h)

  # Check that clusters number is based on k.
  expect_true(all(labels %in% 1:k))
})

test_that("cut_tree throws error if model not fitted", {
  obj <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
  expect_error(obj$cut_tree(k = 2), "Your model is not fitted with any data!")
})

test_that("cut_tree throw and error is no parameter specified", {
  obj <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
  data_quant <- data.frame(
    var1 = rnorm(10),
    var2 = rnorm(10),
    var3 = rnorm(10)
  )
  obj$fit(data_quant)
  expect_error(obj$cut_tree())
})


# ==============================================================================
# Test compute_centroids method
# ==============================================================================

test_that("compute_centroids works with purely quantitative data", {

  # Mock object
  obj <- HClustVar$new(vartype = "quant", dist.metric = "rsquare")
  data_quant <- data.frame(
    var1 = rnorm(10),
    var2 = rnorm(10),
    var3 = rnorm(10)
  )

  obj$fit(data_quant)
  obj$cut_tree(3)

  # Checks
  expect_true(!is.null(obj$centroids))
  expect_true(!is.null(obj$clusters.eigen))
  expect_equal(ncol(obj$centroids), obj$n_clusters)
  expect_true(all(obj$clusters.eigen > 0))
})


test_that("compute_centroids works with purely qualitative data", {
  obj <- HClustVar$new()
  df <- data.frame(
    A = as.factor(sample(letters[1:3], 10, TRUE)),
    B = as.factor(sample(letters[1:2], 10, TRUE)),
    C = as.factor(sample(letters[1:4], 10, TRUE))
  )

  obj$fit(df)
  obj$cut_tree(3)

  expect_true(!is.null(obj$centroids))
  expect_equal(ncol(obj$centroids), obj$n_clusters)
})


test_that("compute_centroids works with mixed data (FAMD case)", {
  obj <- HClustVar$new()
  df <- data.frame(
    A = rnorm(10),
    B = as.factor(sample(letters[1:2], 10, TRUE)),
    C = rnorm(10)
  )


  obj$fit(df)
  obj$cut_tree(3)

  expect_true(!is.null(obj$centroids))
  expect_equal(ncol(obj$centroids), obj$n_clusters)
})


test_that("compute_centroids updates attributes correctly", {
  obj <- HClustVar$new()
  df <- data.frame(
    A = rnorm(10),
    B = as.factor(sample(letters[1:2], 10, TRUE)),
    C = rnorm(10)
  )

  obj$fit(df)
  obj$cut_tree(3)

  expect_true(!is.null(colnames(obj$centroids)))
  expect_equal(length(obj$clusters.eigen), obj$n_clusters)
})


# ==============================================================================
# Test correlation_ratio method
# ==============================================================================

# --- Create a mock instance of HClustVar ---
test_obj <- HClustVar$new(vartype = "quant")

# Access the private method
correlation_ratio <- test_obj$.__enclos_env__$private$correlation_ratio

# --- Sample data for testing ---
quanti <- c(10, 12, 9, 8, 14, 15, 22, 24, 25, 23)
quali <- factor(c("A","A","A","A","B","B","B","B","B","B"))

# --- Unit Tests ---

test_that("correlation_ratio returns a numeric value", {
  eta2 <- correlation_ratio(quali, quanti)
  expect_type(eta2, "double")
})



test_that("correlation_ratio returns near 0 when there is no relationship", {
  # Random data should produce no correlation
  set.seed(1)
  quali_random <- factor(sample(c("A","B","C"), 100, replace = TRUE))
  quanti_random <- rnorm(100)
  eta2 <- correlation_ratio(quali_random, quanti_random)
  expect_true(eta2 < 0.1)
  expect_true(eta2 > 0)
})

test_that("correlation_ratio returns 1 when the separation is perfect", {
  # Perfect group separation should yield eta² = 1
  quanti_perf <- c(1, 2, 3, 4, 10, 11, 12, 13)
  quali_perf <- factor(c("A","A","A","A","B","B","B","B"))
  eta2 <- correlation_ratio(quali_perf, quanti_perf)
  expect_gt(round(eta2, 5), 0.9)
})


# ==============================================================================
# Test predict method
# ==============================================================================

test_that("predict() works for quantitative data", {
  # Create simple quantitative dataset
  df <- data.frame(
    var1 = rnorm(50),
    var2 = rnorm(50),
    var3 = rnorm(50)
  )

  # Fit hierarchical clustering
  model <- HClustVar$new(vartype = "quant", dist.metric = "r")
  model$fit(df)
  model$cut_tree(k = 2)

  # Predict using the same dataset (quantitative)
  preds <- model$predict(df)

  # --- Assertions ---
  expect_type(preds, "integer")          # output should be integer vector
  expect_length(preds, ncol(df))         # one prediction per variable
  expect_true(all(preds %in% 1:model$n_clusters)) # each belongs to a cluster
})


test_that("predict() works for qualitative data", {
  # Create categorical dataset
  df <- data.frame(
    cat1 = factor(sample(letters[1:3], 40, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 40, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 40, replace = TRUE))
  )

  # Fit model for qualitative data
  model <- HClustVar$new(vartype = "qual")
  model$fit(df)
  model$cut_tree(k = 3)

  # Predict with same categorical data
  preds <- model$predict(df)

  # --- Assertions ---
  expect_type(preds, "integer")
  expect_length(preds, ncol(df))
  expect_true(all(preds %in% 1:model$n_clusters))
})


test_that("predict() works for mixed data.", {
  obj <- HClustVar$new()
  df <- data.frame(
    A = rnorm(10),
    B = as.factor(sample(letters[1:2], 10, TRUE)),
    C = rnorm(10)
  )


  obj$fit(df)
  obj$cut_tree(3)

  # predict on the same df.
  preds <- obj$predict(df)

  expect_length(preds, ncol(df))
  expect_true(all(preds %in% 1:obj$n_clusters))
})


test_that("predict() raises error if model is not fitted", {
  model <- HClustVar$new(vartype = "quant")
  df <- data.frame(a = rnorm(10), b = rnorm(10))
  expect_error(model$predict(df), "Model must be fitted")
})


test_that("predict() raises error if tree not cut", {
  model <- HClustVar$new(vartype = "quant")
  df <- data.frame(a = rnorm(10), b = rnorm(10))
  model$fit(df)
  expect_error(model$predict(df), "Tree must be cut")
})


test_that("predict() raises error for incompatible new_data rows", {
  model <- HClustVar$new(vartype = "quant")
  df <- data.frame(a = rnorm(10), b = rnorm(10))
  model$fit(df)
  model$cut_tree(k = 2)

  # new_data with wrong number of rows
  new_df <- data.frame(a = rnorm(5), b = rnorm(5))
  expect_error(model$predict(new_df), "same number of rows")
})
