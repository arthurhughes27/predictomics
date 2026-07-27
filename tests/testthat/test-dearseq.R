# -----------------------------------------------------------------------------
# dearseq selection
# -----------------------------------------------------------------------------

.make_dearseq_data <- function(n_ind = 12, p = 8, seed = 1) {
  set.seed(seed)
  n <- n_ind * 2
  individual_id <- rep(seq_len(n_ind), each = 2)
  timepoint     <- rep(c(0, 1), times = n_ind)
  treatment     <- rep(rep(c(0, 1), each = 2), times = n_ind / 2)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- rnorm(n)
  list(X = X, Y = Y, individual_id = individual_id, timepoint = timepoint,
       treatment = treatment)
}

test_that("selection_params$method = 'dearseq' is a valid option with sensible defaults", {
  d <- .make_dearseq_data()
  params <- list(method = "dearseq", threshold = 0.05,
                 dearseq_mode = "classic")
  expect_error(.validate_selection_params(params, p = ncol(d$X)), NA)
})

test_that("dearseq_mode defaults to 'classic' when not specified", {
  d <- .make_dearseq_data()
  params <- list(method = "dearseq", threshold = 0.05)
  expect_error(.validate_selection_params(params, p = ncol(d$X)), NA)
})

test_that("dearseq-specific params are validated", {
  d <- .make_dearseq_data()
  base <- list(method = "dearseq", threshold = 0.05)

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_mode = "bad")),
                               p = ncol(d$X)),
    "dearseq_mode"
  )

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_which_test = "bad")),
                               p = ncol(d$X)),
    "dearseq_which_test"
  )

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_preprocessed = "yes")),
                               p = ncol(d$X)),
    "dearseq_preprocessed"
  )

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_padjust_methods = "bad")),
                               p = ncol(d$X)),
    "dearseq_padjust_methods"
  )

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_which_weights = "bad")),
                               p = ncol(d$X)),
    "dearseq_which_weights"
  )

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_n_perm = -1)),
                               p = ncol(d$X)),
    "dearseq_n_perm"
  )

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_kernel = "bad")),
                               p = ncol(d$X)),
    "dearseq_kernel"
  )

  expect_error(
    .validate_selection_params(modifyList(base, list(dearseq_bw = TRUE)),
                               p = ncol(d$X)),
    "dearseq_bw"
  )
})

test_that("run_selection errors informatively for dearseq classic mode without treatment", {
  d <- .make_dearseq_data()

  expect_error(
    run_selection(
      X_train = d$X,
      params  = list(method = "dearseq", dearseq_mode = "classic",
                     threshold = 0.05)
    ),
    "treatment must be provided"
  )
})

test_that("run_selection errors when classic mode gets individual_id without timepoint", {
  d <- .make_dearseq_data()

  expect_error(
    run_selection(
      X_train       = d$X,
      treatment     = d$treatment,
      individual_id = d$individual_id,
      params        = list(method = "dearseq", dearseq_mode = "classic",
                           threshold = 0.05)
    ),
    "supplied together"
  )
})

test_that("run_selection errors informatively for dearseq paired mode without pairing", {
  d <- .make_dearseq_data()

  expect_error(
    run_selection(
      X_train = d$X,
      params  = list(method = "dearseq", dearseq_mode = "paired",
                     threshold = 0.05)
    ),
    "individual_id"
  )
})

test_that("run_selection errors on unpaired individual_id/timepoint for dearseq paired mode", {
  d <- .make_dearseq_data()
  timepoint <- d$timepoint
  timepoint[2] <- 0  # individual 1 now has two pre-treatment rows

  expect_error(
    run_selection(
      X_train       = d$X,
      individual_id = d$individual_id,
      timepoint     = timepoint,
      params        = list(method = "dearseq", dearseq_mode = "paired",
                           threshold = 0.05)
    ),
    "one pre-treatment"
  )
})

test_that("run_selection requires the dearseq package once inputs are otherwise valid", {
  testthat::skip_if(requireNamespace("dearseq", quietly = TRUE),
                    "dearseq is installed; skipping missing-package test")
  d <- .make_dearseq_data()

  expect_error(
    run_selection(
      X_train   = d$X,
      treatment = d$treatment,
      params    = list(method = "dearseq", dearseq_mode = "classic",
                       threshold = 0.05)
    ),
    "dearseq"
  )
})

test_that("run_selection computes dearseq scores for all genes (classic mode)", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  res <- run_selection(
    X_train   = d$X,
    treatment = d$treatment,
    params    = list(method = "dearseq", dearseq_mode = "classic",
                     threshold = 1.1)
  )

  expect_equal(length(res$selection_scores), ncol(d$X))
  expect_setequal(names(res$selection_scores), colnames(d$X))
  expect_true(all(diff(res$selection_scores) >= 0))
})

test_that("run_selection computes dearseq scores for all genes (paired mode)", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  res <- run_selection(
    X_train       = d$X,
    individual_id = d$individual_id,
    timepoint     = d$timepoint,
    params        = list(method = "dearseq", dearseq_mode = "paired",
                         threshold = 1.1)
  )

  expect_equal(length(res$selection_scores), ncol(d$X))
  expect_setequal(names(res$selection_scores), colnames(d$X))
})

test_that("predict_cv applies dearseq upfront, before gene_level_fc, regardless of outside_cv", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  result <- predict_cv(
    Y = d$Y, X = d$X,
    folds = 2,
    selection_params    = list(method = "dearseq", dearseq_mode = "classic",
                               threshold = 1.1),
    engineering_params  = list(method = "engineer", gene_level_fc = TRUE),
    individual_id = d$individual_id,
    timepoint     = d$timepoint,
    treatment     = d$treatment,
    outside_cv    = FALSE,
    verbose       = FALSE
  )

  expect_false(is.null(result$dearseq_selection))
  expect_null(result$fold_selection_diagnostics)
  expect_null(result$outside_cv_selection)
  expect_equal(result$n_samples_modelled, length(unique(d$individual_id)))
})
