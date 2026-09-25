# -----------------------------------------------------------------------------
# .compute_relative_gain(): the vectorised (FWL-based) implementation must
# reproduce a naive per-feature lm()-based implementation exactly (up to
# floating-point tolerance) - see R/run_selection.R's Details for the
# Frisch-Waugh-Lovell derivation this equivalence relies on.
# -----------------------------------------------------------------------------

#' Naive, unvectorised reference implementation (one lm() call per feature
#' per inner fold) - a direct transcription of .compute_relative_gain()'s
#' pre-vectorisation logic, kept here only as a correctness oracle for the
#' tests below.
.naive_relative_gain <- function(X_train, Y_train, covariates, metric,
                                  inner_folds, seed) {

  n <- nrow(X_train)
  p <- ncol(X_train)
  feat_names <- colnames(X_train)

  inner_fold_ids <- make_folds(n = n, cv_type = "kfold", k = inner_folds, seed = seed)
  has_covariates <- !is.null(covariates) && ncol(covariates) > 0L

  baseline_pred <- numeric(n)
  for (f in seq_len(inner_folds)) {
    tr <- which(inner_fold_ids != f)
    tst <- which(inner_fold_ids == f)
    Y_tr <- Y_train[tr]
    if (has_covariates) {
      cov_tr <- as.data.frame(covariates[tr, , drop = FALSE])
      cov_tst <- as.data.frame(covariates[tst, , drop = FALSE])
      df_tr <- cbind(data.frame(.Y = Y_tr), cov_tr)
      fit <- lm(.Y ~ ., data = df_tr)
      baseline_pred[tst] <- predict(fit, newdata = cov_tst)
    } else {
      df_tr <- data.frame(.Y = Y_tr)
      fit <- lm(.Y ~ 1, data = df_tr)
      baseline_pred[tst] <- predict(fit, newdata = data.frame(.intercept = rep(1, length(tst))))
    }
  }
  baseline_score <- .compute_metric(Y_train, baseline_pred, metric)

  gains <- numeric(p)
  names(gains) <- feat_names

  for (j in seq_len(p)) {
    feat_pred <- numeric(n)
    for (f in seq_len(inner_folds)) {
      tr <- which(inner_fold_ids != f)
      tst <- which(inner_fold_ids == f)
      Y_tr <- Y_train[tr]
      feat_j <- X_train[, j]
      if (has_covariates) {
        cov_tr <- as.data.frame(covariates[tr, , drop = FALSE])
        cov_tst <- as.data.frame(covariates[tst, , drop = FALSE])
        df_tr <- cbind(data.frame(.Y = Y_tr, .feat = feat_j[tr]), cov_tr)
        df_tst <- cbind(data.frame(.feat = feat_j[tst]), cov_tst)
      } else {
        df_tr <- data.frame(.Y = Y_tr, .feat = feat_j[tr])
        df_tst <- data.frame(.feat = feat_j[tst])
      }
      fit <- lm(.Y ~ ., data = df_tr)
      feat_pred[tst] <- predict(fit, newdata = df_tst)
    }
    feature_score <- .compute_metric(Y_train, feat_pred, metric)
    gains[j] <- .compute_gain(baseline_score, feature_score, metric)
  }

  gains
}

.make_relative_gain_data <- function(n = 40, p = 15, q = 2, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  covariates <- matrix(rnorm(n * q), nrow = n, ncol = q,
                        dimnames = list(NULL, paste0("cov", seq_len(q))))
  Y <- X[, 1] * 2 - X[, 2] + covariates[, 1] * 0.5 + rnorm(n)
  list(X = X, Y = Y, covariates = covariates)
}

test_that(".compute_relative_gain matches the naive per-feature lm() implementation, with covariates", {
  d <- .make_relative_gain_data()

  for (metric in c("rmse", "srmse", "r2", "spearman")) {
    fast  <- .compute_relative_gain(d$X, d$Y, d$covariates, metric = metric,
                                     inner_folds = 5, seed = 123)
    naive <- .naive_relative_gain(d$X, d$Y, d$covariates, metric = metric,
                                   inner_folds = 5, seed = 123)

    expect_equal(fast[names(naive)], naive, tolerance = 1e-8,
                 info = paste("metric =", metric))
  }
})

test_that(".compute_relative_gain matches the naive per-feature lm() implementation, no covariates", {
  d <- .make_relative_gain_data()

  for (metric in c("rmse", "srmse", "r2", "spearman")) {
    fast  <- .compute_relative_gain(d$X, d$Y, covariates = NULL, metric = metric,
                                     inner_folds = 5, seed = 123)
    naive <- .naive_relative_gain(d$X, d$Y, covariates = NULL, metric = metric,
                                   inner_folds = 5, seed = 123)

    expect_equal(fast[names(naive)], naive, tolerance = 1e-8,
                 info = paste("metric =", metric))
  }
})

test_that("run_selection(method = 'relative_gain') still selects sensible top features after vectorisation", {
  d <- .make_relative_gain_data()

  res <- run_selection(
    X_train = d$X, Y_train = d$Y, covariates = d$covariates,
    params = list(method = "relative_gain", top_n = 3,
                  relative_gain_inner_folds = 5, relative_gain_seed = 123)
  )

  # gene1/gene2 are the true signal features (see .make_relative_gain_data());
  # both should rank ahead of at least one pure-noise gene.
  expect_true(all(c("gene1", "gene2") %in% res$selected_features))
})
