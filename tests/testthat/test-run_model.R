# -----------------------------------------------------------------------------
# model_params$impute: NA support in run_model()/predict_model()
# -----------------------------------------------------------------------------

.make_model_data <- function(n = 40, p = 10, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- X[, 1] * 2 + rnorm(n)
  list(X = X, Y = Y)
}

test_that(".impute_matrix computes and applies mean/median fill values", {
  X <- matrix(c(1, 2, NA, 4, NA, 6), nrow = 3, ncol = 2,
             dimnames = list(NULL, c("a", "b")))

  res_mean <- .impute_matrix(X, method = "mean")
  expect_equal(res_mean$values[["a"]], mean(c(1, 2), na.rm = TRUE))
  expect_equal(res_mean$values[["b"]], mean(c(4, 6), na.rm = TRUE))
  expect_false(anyNA(res_mean$X))

  res_median <- .impute_matrix(X, method = "median")
  expect_equal(res_median$values[["a"]], median(c(1, 2)))
  expect_false(anyNA(res_median$X))
})

test_that(".impute_matrix reapplies supplied fill values without recomputing", {
  X_train <- matrix(c(1, 2, NA, 10, 20, 30), nrow = 3, ncol = 2,
                    dimnames = list(NULL, c("a", "b")))
  train_res <- .impute_matrix(X_train, method = "mean")

  X_test <- matrix(c(NA, 5), nrow = 1, ncol = 2,
                   dimnames = list(NULL, c("a", "b")))
  test_res <- .impute_matrix(X_test, method = "mean",
                             fill_values = train_res$values)

  expect_equal(test_res$X[1, "a"], train_res$values[["a"]])
  expect_equal(test_res$X[1, "b"], 5)
})

test_that(".impute_matrix errors informatively on a fully-NA column", {
  X <- matrix(c(NA, NA, 1, 2), nrow = 2, ncol = 2,
             dimnames = list(NULL, c("a", "b")))
  expect_error(.impute_matrix(X, method = "mean"), "all values are NA")
})

test_that("run_model errors informatively when X_train has NA and impute = 'none'", {
  d <- .make_model_data()
  X_na <- d$X
  X_na[1, 1] <- NA

  expect_error(
    run_model(X_train = X_na, Y_train = d$Y, params = list(method = "lm")),
    "impute"
  )
})

test_that("run_model + predict_model support NA in X via impute = 'mean'", {
  d <- .make_model_data()
  X_train_na <- d$X
  X_train_na[1, 1] <- NA

  fit <- run_model(X_train = X_train_na, Y_train = d$Y,
                   params = list(method = "lm", impute = "mean"))

  expect_equal(fit$impute_method, "mean")
  expect_false(is.null(fit$impute_values))

  X_test_na <- d$X[1:5, , drop = FALSE]
  X_test_na[2, 1] <- NA

  preds <- predict_model(fit, X_new = X_test_na)
  expect_equal(length(preds), 5L)
  expect_false(anyNA(preds))
})

test_that("run_model + predict_model support NA in X via impute = 'median' for glmnet", {
  d <- .make_model_data()
  X_train_na <- d$X
  X_train_na[1, 1] <- NA
  X_train_na[2, 3] <- NA

  fit <- run_model(X_train = X_train_na, Y_train = d$Y,
                   params = list(method = "glmnet", inner_folds = 3,
                                 impute = "median"))

  X_test_na <- d$X[1:5, , drop = FALSE]
  X_test_na[1, 2] <- NA

  preds <- predict_model(fit, X_new = X_test_na)
  expect_equal(length(preds), 5L)
  expect_false(anyNA(preds))
})

test_that("predict_model errors informatively when X_new has NA but impute = 'none'", {
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y, params = list(method = "lm"))

  X_test_na <- d$X[1:5, , drop = FALSE]
  X_test_na[1, 1] <- NA

  expect_error(predict_model(fit, X_new = X_test_na), "impute")
})

test_that("predict_cv end-to-end supports NA in X via model_params$impute", {
  d <- .make_model_data(n = 40, p = 10)
  X_na <- d$X
  set.seed(2)
  na_idx <- cbind(sample(nrow(X_na), 5), sample(ncol(X_na), 5))
  X_na[na_idx] <- NA

  result <- predict_cv(
    Y = d$Y, X = X_na,
    model_params = list(method = "lm", impute = "mean"),
    folds   = 4,
    verbose = FALSE
  )

  expect_equal(length(result$predicted), nrow(d$X))
  expect_false(anyNA(result$predicted))
})

test_that("predict_cv errors informatively when X has NA and impute = 'none' (default)", {
  d <- .make_model_data(n = 40, p = 10)
  X_na <- d$X
  X_na[1, 1] <- NA

  expect_error(
    predict_cv(Y = d$Y, X = X_na, folds = 4, verbose = FALSE),
    "impute"
  )
})

# -----------------------------------------------------------------------------
# model_params$scale
# -----------------------------------------------------------------------------

test_that(".validate_model_params rejects a non-logical/NA scale", {
  expect_error(
    .validate_model_params(list(method = "lm", scale = "yes")),
    "model_params\\$scale"
  )
  expect_error(
    .validate_model_params(list(method = "lm", scale = c(TRUE, FALSE))),
    "model_params\\$scale"
  )
  expect_error(
    .validate_model_params(list(method = "lm", scale = NA)),
    "model_params\\$scale"
  )
})

test_that(".validate_model_params accepts scale = TRUE/FALSE and defaults it", {
  expect_silent(.validate_model_params(list(method = "lm", scale = TRUE)))
  expect_silent(.validate_model_params(list(method = "lm", scale = FALSE)))
  expect_silent(.validate_model_params(list(method = "lm")))
})

test_that("run_model records scale = FALSE by default", {
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y, params = list(method = "lm"))
  expect_identical(fit$scale, FALSE)
})

test_that("run_model + predict_model give near-identical lm predictions with scale = TRUE/FALSE", {
  d <- .make_model_data()

  fit_unscaled <- run_model(X_train = d$X, Y_train = d$Y,
                            params = list(method = "lm", scale = FALSE))
  fit_scaled <- run_model(X_train = d$X, Y_train = d$Y,
                          params = list(method = "lm", scale = TRUE))

  expect_identical(fit_scaled$scale, TRUE)

  X_test <- d$X[1:5, , drop = FALSE]
  preds_unscaled <- predict_model(fit_unscaled, X_new = X_test)
  preds_scaled <- predict_model(fit_scaled, X_new = X_test)

  # OLS predictions are invariant to affine per-feature rescaling of the
  # predictors, so the two should agree up to numerical tolerance.
  expect_equal(preds_scaled, preds_unscaled, tolerance = 1e-6)
})

test_that("run_model + predict_model support scale = TRUE for glmnet", {
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y,
                   params = list(method = "glmnet", inner_folds = 3,
                                 scale = TRUE))
  expect_identical(fit$scale, TRUE)

  preds <- predict_model(fit, X_new = d$X[1:5, , drop = FALSE])
  expect_equal(length(preds), 5L)
  expect_false(anyNA(preds))
})

test_that("run_model + predict_model support scale = TRUE for ranger", {
  testthat::skip_if_not_installed("ranger")
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y,
                   params = list(method = "ranger", inner_folds = 3,
                                 scale = TRUE))
  expect_identical(fit$scale, TRUE)

  preds <- predict_model(fit, X_new = d$X[1:5, , drop = FALSE])
  expect_equal(length(preds), 5L)
  expect_false(anyNA(preds))
})

test_that("run_model + predict_model support scale = TRUE for svr", {
  testthat::skip_if_not_installed("kernlab")
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y,
                   params = list(method = "svr", inner_folds = 3,
                                 scale = TRUE))
  expect_identical(fit$scale, TRUE)

  preds <- predict_model(fit, X_new = d$X[1:5, , drop = FALSE])
  expect_equal(length(preds), 5L)
  expect_false(anyNA(preds))
})

# -----------------------------------------------------------------------------
# model_params$compute_importance
# -----------------------------------------------------------------------------

test_that(".validate_model_params rejects a non-logical/NA compute_importance", {
  expect_error(
    .validate_model_params(list(method = "lm", compute_importance = "yes")),
    "model_params\\$compute_importance"
  )
  expect_error(
    .validate_model_params(list(method = "lm",
                                compute_importance = c(TRUE, FALSE))),
    "model_params\\$compute_importance"
  )
  expect_error(
    .validate_model_params(list(method = "lm", compute_importance = NA)),
    "model_params\\$compute_importance"
  )
})

test_that("run_model does not compute feature_importance by default", {
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y, params = list(method = "lm"))
  expect_null(fit$feature_importance)
  expect_true(is.na(fit$importance_type))
})

test_that("run_model computes feature_importance for lm covering all features", {
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y,
                   params = list(method = "lm", scale = TRUE,
                                 compute_importance = TRUE))

  expect_identical(fit$importance_type, "coefficient")
  expect_setequal(names(fit$feature_importance), colnames(d$X))
  # sorted by decreasing absolute value
  expect_equal(fit$feature_importance,
              fit$feature_importance[order(abs(fit$feature_importance),
                                           decreasing = TRUE)])
})

test_that("run_model computes feature_importance for glmnet (dense, including ridge)", {
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y,
                   params = list(method = "ridge", inner_folds = 3,
                                 scale = TRUE, compute_importance = TRUE))

  expect_identical(fit$importance_type, "coefficient")
  expect_setequal(names(fit$feature_importance), colnames(d$X))
})

test_that("run_model computes non-negative feature_importance for ranger", {
  testthat::skip_if_not_installed("ranger")
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y,
                   params = list(method = "ranger", inner_folds = 3,
                                 compute_importance = TRUE))

  expect_identical(fit$importance_type, "impurity")
  expect_setequal(names(fit$feature_importance), colnames(d$X))
  expect_true(all(fit$feature_importance >= 0))
})

test_that("run_model computes feature_importance for svr via the linear weight vector", {
  testthat::skip_if_not_installed("kernlab")
  d <- .make_model_data()
  fit <- run_model(X_train = d$X, Y_train = d$Y,
                   params = list(method = "svr", inner_folds = 3,
                                 scale = TRUE, compute_importance = TRUE))

  expect_identical(fit$importance_type, "coefficient")
  expect_setequal(names(fit$feature_importance), colnames(d$X))
  # Regression test: kernlab's own xmatrix slot does not reliably carry
  # column names, so names must come from positional mapping, never from
  # kernlab's self-reported (possibly NA/missing) names.
  expect_false(anyNA(names(fit$feature_importance)))
  expect_equal(length(fit$feature_importance), ncol(d$X))
})

test_that("predict_cv warns when compute_importance = TRUE, coefficient-based method, scale != TRUE", {
  d <- .make_model_data()
  expect_warning(
    predict_cv(Y = d$Y, X = d$X,
              model_params = list(method = "lm", compute_importance = TRUE),
              folds = 3, verbose = FALSE),
    "model_params\\$scale"
  )
})

test_that("predict_cv does not warn about scale when compute_importance + scale = TRUE", {
  d <- .make_model_data()
  warnings_seen <- character(0)
  withCallingHandlers(
    predict_cv(Y = d$Y, X = d$X,
              model_params = list(method = "lm", scale = TRUE,
                                  compute_importance = TRUE),
              folds = 3, verbose = FALSE),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("model_params\\$scale", warnings_seen)))
})

test_that("predict_cv populates fold_feature_importance only when requested", {
  d <- .make_model_data()

  result_off <- predict_cv(Y = d$Y, X = d$X,
                           model_params = list(method = "lm"),
                           folds = 3, verbose = FALSE)
  expect_null(result_off$fold_feature_importance)

  result_on <- suppressWarnings(predict_cv(
    Y = d$Y, X = d$X,
    model_params = list(method = "lm", compute_importance = TRUE),
    folds = 3, verbose = FALSE
  ))
  expect_equal(length(result_on$fold_feature_importance), 3L)
  expect_setequal(names(result_on$fold_feature_importance[[1]]$scores),
                  colnames(d$X))
})
