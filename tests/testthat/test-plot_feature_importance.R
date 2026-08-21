# -----------------------------------------------------------------------------
# plot_feature_importance()
# -----------------------------------------------------------------------------

.make_importance_data <- function(n = 30, p = 10, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- X[, 1] * 2 + rnorm(n)
  list(X = X, Y = Y)
}

test_that("plot_feature_importance errors without a predictomics object", {
  expect_error(plot_feature_importance(list(a = 1)), "predictomics object")
})

test_that("plot_feature_importance errors when compute_importance was not requested", {
  d <- .make_importance_data()
  result <- predict_cv(Y = d$Y, X = d$X, model_params = list(method = "lm"),
                       folds = 3, verbose = FALSE)

  expect_error(plot_feature_importance(result), "compute_importance")
})

test_that("plot_feature_importance errors on a non-positive top_n", {
  d <- .make_importance_data()
  result <- suppressWarnings(predict_cv(
    Y = d$Y, X = d$X,
    model_params = list(method = "lm", compute_importance = TRUE),
    folds = 3, verbose = FALSE
  ))

  expect_error(plot_feature_importance(result, top_n = 0), "top_n")
})

test_that("plot_feature_importance produces a ggplot ranked by mean |importance|", {
  d <- .make_importance_data()
  result <- predict_cv(
    Y = d$Y, X = d$X,
    model_params = list(method = "lm", scale = TRUE, compute_importance = TRUE),
    folds = 3, verbose = FALSE
  )

  p <- plot_feature_importance(result, top_n = 5)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 5L)
  expect_equal(p$data$importance, sort(p$data$importance, decreasing = TRUE))
})

test_that("plot_feature_importance warns when scale != TRUE for a coefficient-based method", {
  d <- .make_importance_data()
  result <- suppressWarnings(predict_cv(
    Y = d$Y, X = d$X,
    model_params = list(method = "lm", compute_importance = TRUE),
    folds = 3, verbose = FALSE
  ))

  expect_warning(plot_feature_importance(result), "scale")
})

test_that("plot_feature_importance does not warn for ranger (impurity-based)", {
  testthat::skip_if_not_installed("ranger")
  d <- .make_importance_data()
  result <- predict_cv(
    Y = d$Y, X = d$X,
    model_params = list(method = "ranger", inner_folds = 3,
                        compute_importance = TRUE),
    folds = 3, verbose = FALSE
  )

  warnings_seen <- character(0)
  withCallingHandlers(
    plot_feature_importance(result),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("scale", warnings_seen)))
})
