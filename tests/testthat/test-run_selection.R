# -----------------------------------------------------------------------------
# run_selection(): top_n exceeding the number of available features
# -----------------------------------------------------------------------------

.make_selection_data <- function(n = 20, p = 6, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- X[, 1] * 2 + rnorm(n)
  list(X = X, Y = Y)
}

test_that(".validate_selection_params no longer errors when top_n exceeds the feature count", {
  d <- .make_selection_data()

  expect_error(
    .validate_selection_params(
      list(method = "variance", top_n = 5000), p = ncol(d$X)
    ),
    NA
  )
})

test_that("run_selection selects all features (with a message) when top_n exceeds the feature count", {
  d <- .make_selection_data()

  msgs <- testthat::capture_messages(
    res <- run_selection(
      X_train = d$X,
      params  = list(method = "variance", top_n = 5000)
    )
  )

  expect_true(any(grepl("exceeds \\(or equals\\)", msgs)))
  expect_true(any(grepl("without computing selection scores", msgs)))
  expect_setequal(res$selected_features, colnames(d$X))
  expect_true(all(is.na(res$selection_scores)))
  expect_equal(res$top_n, ncol(d$X))
})

test_that("run_selection still computes real scores when top_n is within range", {
  d <- .make_selection_data()

  msgs <- testthat::capture_messages(
    res <- run_selection(
      X_train = d$X, Y_train = d$Y,
      params  = list(method = "pearson", top_n = 3)
    )
  )

  expect_false(any(grepl("exceeds \\(or equals\\)", msgs)))
  expect_equal(length(res$selected_features), 3L)
  expect_false(any(is.na(res$selection_scores)))
})

test_that("run_selection selects all features when top_n exactly equals the feature count", {
  d <- .make_selection_data()

  res <- run_selection(
    X_train = d$X,
    params  = list(method = "variance", top_n = ncol(d$X))
  )

  expect_setequal(res$selected_features, colnames(d$X))
})

test_that("relative_gain's negative-gain check is skipped when top_n exceeds the feature count", {
  d <- .make_selection_data()

  expect_error(
    run_selection(
      X_train = d$X, Y_train = d$Y,
      params  = list(method = "relative_gain", top_n = 5000,
                     relative_gain_inner_folds = 2)
    ),
    NA
  )
})

test_that("predict_cv completes when selection top_n exceeds the number of available features", {
  d <- .make_selection_data()

  result <- predict_cv(
    Y = d$Y, X = d$X,
    selection_params = list(method = "variance", top_n = 5000),
    folds   = 4,
    verbose = FALSE
  )

  expect_equal(length(result$observed), nrow(d$X))
})
