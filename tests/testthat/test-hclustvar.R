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
    regexp = "dist.metric ignored for qualitative data"
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

# One column dataframe
df02 <- df01[1]

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

test_that("quantile_discretisation works on single column dataframe", {
  df_disc <- quantile_discretisation(df02, quanti_index = 1, n_groups = 4)

  # Quantitatives should become factors.
  expect_true(is.factor(df_disc[, 1]))

  # We should have Q1 to Q4 labels.
  expect_equal(levels(df_disc[[1]]), c("Q1", "Q2", "Q3", "Q4"))
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
  expect_error(obj$cut_tree(k = 2), "Your model should be fitted on data first.")
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


# ============================================================================
# Predict method
# ============================================================================

library(testthat)

# ============================================================================
# 1. Prerequisites validation
# ============================================================================

test_that("predict raise an error if not fitted", {
  hc <- HClustVar$new(vartype = "quant")
  new_data <- data.frame(new_var = rnorm(50))

  expect_error(
    hc$predict(new_data),
    "Model must be fitted before prediction"
  )
})

test_that("Predict raise an error if the tree not cut.", {
  set.seed(123)
  data <- data.frame(v1 = rnorm(50), v2 = rnorm(50), v3 = rnorm(50))
  new_data <- data.frame(new_var = rnorm(50))

  hc <- HClustVar$new(vartype = "quant")
  hc$fit(data)

  expect_error(
    hc$predict(new_data),
    "Tree must be cut before prediction"
  )
})

test_that("predict raise an error on non df or matrix.", {
  set.seed(123)
  data <- data.frame(v1 = rnorm(50), v2 = rnorm(50), v3 = rnorm(50))

  hc <- HClustVar$new(vartype = "quant")
  hc$fit(data)
  hc$cut_tree(k = 2)

  expect_error(
    hc$predict(list(a = 1, b = 2)),
    "new_data must be a data.frame or matrix"
  )

  expect_error(
    hc$predict(c(1, 2, 3)),
    "new_data must be a data.frame or matrix"
  )
})

test_that("predict raise an error if different number of lines", {
  set.seed(123)
  data <- data.frame(v1 = rnorm(50), v2 = rnorm(50), v3 = rnorm(50))
  wrong_data <- data.frame(new_var = rnorm(100))

  hc <- HClustVar$new(vartype = "quant")
  hc$fit(data)
  hc$cut_tree(k = 2)

  expect_error(
    hc$predict(wrong_data),
    "new_data must have the same number of rows as training data"
  )
})

test_that("predict raise an error if prediction on empty data", {
  set.seed(123)
  data <- data.frame(v1 = rnorm(50), v2 = rnorm(50), v3 = rnorm(50))
  empty_data <- data.frame()

  hc <- HClustVar$new(vartype = "quant")
  hc$fit(data)
  hc$cut_tree(k = 2)

  expect_error(
    hc$predict(empty_data)
  )
})

test_that("predict work well on matrices.", {
  set.seed(123)
  data <- data.frame(v1 = rnorm(50), v2 = rnorm(50), v3 = rnorm(50))
  new_data_matrix <- matrix(rnorm(50), ncol = 1)
  colnames(new_data_matrix) <- "new_v"

  hc <- HClustVar$new(vartype = "quant")
  hc$fit(data)
  hc$cut_tree(k = 2)

  expect_no_error(result <- hc$predict(new_data_matrix))
  expect_type(result, "integer")
  expect_named(result, "new_v")
})

# ============================================================================
# 2. Quantitatives variables - WARD method
# ============================================================================

test_that("predict - quantitatives Variables (ward.D, metric = rsquare)", {
  set.seed(456)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100) + 2,
    v3 = rnorm(100),
    v4 = rnorm(100) + 2
  )

  new_data <- data.frame(
    new_v1 = rnorm(100),
    new_v2 = rnorm(100) + 2
  )

  hc <- HClustVar$new(vartype = "quant", dist.metric = "rsquare", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_length(result, ncol(new_data))
  expect_true(all(result %in% 1:2))
  expect_named(result, colnames(new_data))
})

test_that("predict - quantitatives Variables (ward.D, metric = r)", {
  set.seed(789)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", dist.metric = "r", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_length(result, 1)
  expect_true(result %in% 1:2)
})

test_that("predict - quantitatives Variables (ward.D2)", {
  set.seed(111)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100),
    v4 = rnorm(100)
  )

  new_data <- data.frame(
    new_v1 = rnorm(100),
    new_v2 = rnorm(100)
  )

  hc <- HClustVar$new(vartype = "quant", cah.method = "ward.D2")
  hc$fit(data)
  hc$cut_tree(k = 3)

  result <- hc$predict(new_data)

  expect_length(result, 2)
  expect_true(all(result %in% 1:3))
})

test_that("predict - quantitatives Variables (centroid)", {
  set.seed(222)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", cah.method = "centroid")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

# ============================================================================
# 3. Qualitative variables - WARD METHOD
# ============================================================================

test_that("predict - Qualitatives Variables (ward.D)", {
  set.seed(333)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "qual", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_length(result, 1)
  expect_true(result %in% 1:2)
  expect_named(result, "new_cat")
})

test_that("predict - Qualitatives Variables with qualitative variable (ward.D2)", {
  set.seed(444)
  data <- data.frame(
    cat1 = factor(sample(letters[1:4], 100, replace = TRUE)),
    cat2 = factor(sample(letters[5:8], 100, replace = TRUE)),
    cat3 = factor(sample(letters[9:12], 100, replace = TRUE)),
    cat4 = factor(sample(letters[13:16], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    new_cat2 = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "qual", cah.method = "ward.D2")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_length(result, 2)
  expect_true(all(result %in% 1:2))
})

# ============================================================================
# 4. Tests for mixed variables - Ward method
# ============================================================================

test_that("predict - mixed quantitative variables (ward.D)", {
  set.seed(555)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  new_data <- data.frame(new_quant = rnorm(100))

  hc <- HClustVar$new(vartype = "mixed", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_length(result, 1)
  expect_true(result %in% 1:2)
})

test_that("predict - mixed quantitative variables (ward.D)", {
  set.seed(666)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "mixed", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - mixed quantitative variables (ward.D)", {
  set.seed(777)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_quant = rnorm(100),
    new_cat = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "mixed", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_length(result, 2)
  expect_true(all(result %in% 1:2))
  expect_named(result, c("new_quant", "new_cat"))
})

# ============================================================================
# 5. TESTS CORRELATION RATIO
# ============================================================================

test_that("predict - Quantitatives variables with qualitatives data (ward.D)", {
  set.seed(888)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  new_data <- data.frame(new_quant = rnorm(100))

  hc <- HClustVar$new(vartype = "qual", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_length(result, 1)
  expect_true(result %in% 1:2)
})

test_that("predict - Qualitative variable with quantitative data (ward.D)", {
  set.seed(999)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_length(result, 1)
  expect_true(result %in% 1:2)
})

# ============================================================================
# 6. OTHER CAH METHODS - QUANTITATIVES
# ============================================================================

test_that("predict - quantitatives with single method", {
  set.seed(1001)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", dist.metric = "r", cah.method = "single")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_length(result, 1)
  expect_true(result %in% 1:2)
})

test_that("predict - complete method with quantitatives", {
  set.seed(1002)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", dist.metric = "rsquare", cah.method = "complete")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - median method with quantitatives", {
  set.seed(1003)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100),
    v4 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", cah.method = "median")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - average method with quantitatives", {
  set.seed(1004)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", cah.method = "average")
  hc$fit(data)
  hc$cut_tree(k = 3)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:3)
})

test_that("predict - McQuitty method with quantitatives", {
  set.seed(1005)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100),
    v4 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", cah.method = "mcquitty")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

# ============================================================================
# 7. OTHER CAH METHODS - QUALITATIVES
# ============================================================================

test_that("predict - Single method - qualitatives", {
  set.seed(1006)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "qual", cah.method = "single")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - Complete method - qualitatives", {
  set.seed(1007)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "qual", cah.method = "complete")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - Average method - qualitatives", {
  set.seed(1008)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[10:12], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "qual", cah.method = "average")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

# ============================================================================
# 8. OTHER CAH METHODS - MIXED VARIABLES
# ============================================================================

test_that("predict - Single method - mixed variables", {
  set.seed(1009)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_quant = rnorm(100),
    new_cat = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "mixed", cah.method = "single")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_length(result, 2)
  expect_true(all(result %in% 1:2))
})

test_that("predict - Complete method - mixed variables", {
  set.seed(1010)
  data <- data.frame(
    v1 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  new_data <- data.frame(new_quant = rnorm(100))

  hc <- HClustVar$new(vartype = "mixed", cah.method = "complete")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

# ============================================================================
# 9. Correlation ratio with other meethods
# ============================================================================

test_that("predict - quantitative new data with qualitative data (average)", {
  set.seed(1011)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  new_data <- data.frame(new_quant = rnorm(100))

  hc <- HClustVar$new(vartype = "qual", cah.method = "average")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - qualitative new data with quantitative data (average)", {
  set.seed(1012)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "quant", cah.method = "single")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

# ============================================================================
# 10. TESTS WITH MULTIPLE VARIABLES
# ============================================================================

test_that("predict - quantitative variables", {
  set.seed(1013)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100),
    v4 = rnorm(100)
  )

  new_data <- data.frame(
    new_v1 = rnorm(100),
    new_v2 = rnorm(100),
    new_v3 = rnorm(100)
  )

  hc <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 3)

  result <- hc$predict(new_data)

  expect_length(result, 3)
  expect_named(result, c("new_v1", "new_v2", "new_v3"))
  expect_true(all(result %in% 1:3))
})

test_that("predict - qualitative variables", {
  set.seed(1014)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat1 = factor(sample(letters[10:12], 100, replace = TRUE)),
    new_cat2 = factor(sample(letters[13:15], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "qual", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_length(result, 2)
  expect_named(result, c("new_cat1", "new_cat2"))
  expect_true(all(result %in% 1:2))
})

test_that("predict - mixed variables", {
  set.seed(1015)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_v1 = rnorm(100),
    new_cat1 = factor(sample(letters[7:9], 100, replace = TRUE)),
    new_v2 = rnorm(100)
  )

  hc <- HClustVar$new(vartype = "mixed", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  result <- hc$predict(new_data)

  expect_length(result, 3)
  expect_named(result, c("new_v1", "new_cat1", "new_v2"))
  expect_true(all(result %in% 1:2))
})

# ============================================================================
# 11. TESTS PREDICTIONS CONSISTENCY
# ============================================================================

test_that("predict - Consistency with multiple clusters", {
  set.seed(1017)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100),
    v4 = rnorm(100),
    v5 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  # Test avec k = 2
  hc1 <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc1$fit(data)
  hc1$cut_tree(k = 2)
  result1 <- hc1$predict(new_data)
  expect_true(result1 %in% 1:2)

  # Test avec k = 3
  hc2 <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc2$fit(data)
  hc2$cut_tree(k = 3)
  result2 <- hc2$predict(new_data)
  expect_true(result2 %in% 1:3)

  # Test avec k = 4
  hc3 <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc3$fit(data)
  hc3$cut_tree(k = 4)
  result3 <- hc3$predict(new_data)
  expect_true(result3 %in% 1:4)
})

# ============================================================================
# 12. VARTYPE = AUTO tests
# ============================================================================

test_that("predict - vartype=auto detect quantitatives", {
  set.seed(1018)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "auto", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  expect_equal(hc$vartype, "quant")

  result <- hc$predict(new_data)
  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - vartype=auto detect qualitatives", {
  set.seed(1019)
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "auto", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  expect_equal(hc$vartype, "qual")

  result <- hc$predict(new_data)
  expect_type(result, "integer")
  expect_true(result %in% 1:2)
})

test_that("predict - vartype=auto detect mixed data", {
  set.seed(1020)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_quant = rnorm(100),
    new_cat = factor(sample(letters[4:6], 100, replace = TRUE))
  )

  hc <- HClustVar$new(vartype = "auto", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  expect_equal(hc$vartype, "mixed")

  result <- hc$predict(new_data)
  expect_length(result, 2)
  expect_true(all(result %in% 1:2))
})

# ============================================================================
# 13. CORNERCASES
# ============================================================================

test_that("predict - work with only one cluster", {
  set.seed(1021)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 1)

  result <- hc$predict(new_data)

  expect_equal(result[[1]], 1)
})

test_that("predict - work with variables = clusters", {
  set.seed(1022)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data <- data.frame(new_v = rnorm(100))

  hc <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 3)

  result <- hc$predict(new_data)

  expect_true(result %in% 1:3)
})

test_that("predict - handle missing values", {
  set.seed(1023)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  new_data_with_na <- data.frame(new_v = c(rnorm(90), rep(NA, 10)))

  hc <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  # Should work because : use = "complete.obs" in cor()
  expect_no_error(result <- hc$predict(new_data_with_na))
  expect_type(result, "integer")
})

# ============================================================================
# 14. STABILITY TEST
# ============================================================================

test_that("predict - Reproducible results withy same seed", {
  set.seed(1024)
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    v3 = rnorm(100)
  )

  hc <- HClustVar$new(vartype = "quant", cah.method = "ward.D")
  hc$fit(data)
  hc$cut_tree(k = 2)

  set.seed(2024)
  new_data1 <- data.frame(new_v = rnorm(100))
  result1 <- hc$predict(new_data1)

  set.seed(2024)
  new_data2 <- data.frame(new_v = rnorm(100))
  result2 <- hc$predict(new_data2)

  expect_equal(result1, result2)
})


# ============================================================================
# 15. INTEGRATION TESTS
# ============================================================================

test_that("predict - Complete qualitative workflow", {
  set.seed(1026)

  # Création des données
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100) + 2,
    v3 = rnorm(100),
    v4 = rnorm(100) + 2,
    v5 = rnorm(100)
  )

  new_data <- data.frame(
    new_v1 = rnorm(100),
    new_v2 = rnorm(100) + 2
  )


  hc <- HClustVar$new(vartype = "quant", dist.metric = "rsquare", cah.method = "ward.D")
  hc$fit(data)
  labels <- hc$cut_tree(k = 3)
  predictions <- hc$predict(new_data)

  # checks
  expect_true(hc$fitted)
  expect_equal(hc$n_clusters, 3)
  expect_length(labels, 5)
  expect_length(predictions, 2)
  expect_true(all(predictions %in% 1:3))
  expect_named(predictions, c("new_v1", "new_v2"))
})

test_that("predict -  Complete quantitative workflow", {
  set.seed(1027)

  # Data creation
  data <- data.frame(
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    cat3 = factor(sample(letters[7:9], 100, replace = TRUE)),
    cat4 = factor(sample(letters[10:12], 100, replace = TRUE))
  )

  new_data <- data.frame(
    new_cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    new_cat2 = factor(sample(letters[13:15], 100, replace = TRUE))
  )

  # complete workflow
  hc <- HClustVar$new(vartype = "qual", cah.method = "ward.D")
  hc$fit(data)
  labels <- hc$cut_tree(k = 2)
  predictions <- hc$predict(new_data)

  # checks
  expect_true(hc$fitted)
  expect_equal(hc$n_clusters, 2)
  expect_length(labels, 4)
  expect_length(predictions, 2)
  expect_true(all(predictions %in% 1:2))
})

test_that("predict -  Complete mixed workflow", {
  set.seed(1028)

  # data creation
  data <- data.frame(
    v1 = rnorm(100),
    v2 = rnorm(100),
    cat1 = factor(sample(letters[1:3], 100, replace = TRUE)),
    cat2 = factor(sample(letters[4:6], 100, replace = TRUE)),
    v3 = rnorm(100)
  )

  new_data <- data.frame(
    new_v = rnorm(100),
    new_cat = factor(sample(letters[7:9], 100, replace = TRUE))
  )

  # Complete workflow
  hc <- HClustVar$new(vartype = "mixed", cah.method = "ward.D")
  hc$fit(data)
  labels <- hc$cut_tree(k = 3)
  predictions <- hc$predict(new_data)

  # Checks
  expect_true(hc$fitted)
  expect_equal(hc$n_clusters, 3)
  expect_length(labels, 5)
  expect_length(predictions, 2)
  expect_true(all(predictions %in% 1:3))
  expect_named(predictions, c("new_v", "new_cat"))
})
