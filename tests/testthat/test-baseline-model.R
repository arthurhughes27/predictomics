# -----------------------------------------------------------------------------
# Baseline model (X = NULL)
# -----------------------------------------------------------------------------

test_that("X = NULL with no covariates predicts the training-fold mean of Y", {
  set.seed(1)
  n <- 20
  Y <- rnorm(n)

  result <- predict_cv(Y = Y, X = NULL, folds = 4, verbose = FALSE)

  expect_true(result$baseline_model)
  expect_equal(result$n_features_input, 0L)
  expect_equal(length(result$predicted), n)

  fold_ids <- result$fold_ids
  for (k in unique(fold_ids)) {
    train_mean <- mean(Y[fold_ids != k])
    expect_equal(unique(result$predicted[fold_ids == k]), train_mean)
  }
})

test_that("X = NULL with covariates predicts using only covariates", {
  set.seed(1)
  n <- 30
  age <- rnorm(n)
  Y   <- 2 * age + rnorm(n, sd = 0.1)
  covariates <- data.frame(age = age)

  result <- predict_cv(Y = Y, X = NULL, covariates = covariates,
                       folds = 5, verbose = FALSE)

  expect_true(result$baseline_model)
  expect_equal(result$n_features_input, 0L)
  # Predictions should track Y reasonably well since age is highly predictive
  expect_gt(cor(result$observed, result$predicted), 0.8)
})

test_that("X = NULL with treatment_predictor = TRUE uses treatment as a predictor", {
  set.seed(1)
  n <- 30
  treatment <- rep(c(0, 1), each = n / 2)
  Y <- treatment * 5 + rnorm(n, sd = 0.1)

  result <- predict_cv(Y = Y, X = NULL, treatment = treatment,
                       treatment_predictor = TRUE, folds = 5, verbose = FALSE)

  expect_true(result$baseline_model)
  expect_gt(cor(result$observed, result$predicted), 0.8)
})

test_that("X = NULL errors informatively when combined with engineering_params or selection_params", {
  set.seed(1)
  Y <- rnorm(20)

  expect_error(
    predict_cv(Y = Y, X = NULL,
              engineering_params = list(method = "engineer", col_transform = "z"),
              verbose = FALSE),
    "not compatible"
  )

  expect_error(
    predict_cv(Y = Y, X = NULL,
              selection_params = list(method = "variance", top_n = 5),
              verbose = FALSE),
    "not compatible"
  )
})

test_that("X = NULL prints an informative message about the baseline model", {
  set.seed(1)
  Y <- rnorm(20)

  msgs_no_cov <- testthat::capture_messages(
    predict_cv(Y = Y, X = NULL, verbose = TRUE)
  )
  expect_true(any(grepl("baseline model", msgs_no_cov)))

  covariates <- data.frame(age = rnorm(20))
  msgs_with_cov <- testthat::capture_messages(
    predict_cv(Y = Y, X = NULL, covariates = covariates, verbose = TRUE)
  )
  expect_true(any(grepl("baseline model", msgs_with_cov)))
})

test_that("X = NULL baseline model works with cv_type = 'loo'", {
  set.seed(1)
  n <- 12
  Y <- rnorm(n)

  result <- predict_cv(Y = Y, X = NULL, cv_type = "loo", verbose = FALSE)

  expect_true(result$baseline_model)
  expect_equal(result$n_folds, n)
})
