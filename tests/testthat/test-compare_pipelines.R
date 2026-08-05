# -----------------------------------------------------------------------------
# compare_pipelines()
# -----------------------------------------------------------------------------

.make_compare_data <- function(n = 40, p = 20, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- X[, 1] * 2 + rnorm(n)
  list(X = X, Y = Y)
}

test_that("compare_pipelines runs a basic selection comparison", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "selection",
    option_choices = list(
      list(method = "spearman", top_n = 5),
      list(method = "variance", top_n = 5)
    ),
    reference_params = list(
      selection_params = list(method = "pearson", top_n = 5),
      model_params      = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  expect_s3_class(cmp, "predictomics_comparison")
  expect_equal(nrow(cmp$results), 4L)
  expect_setequal(cmp$results$role, c("baseline", "reference", "option"))
  expect_setequal(cmp$results$pipeline,
                  c("Baseline", "Reference", "spearman", "variance"))
  expect_true(all(c("RMSE", "sRMSE", "R2", "SpearmanR") %in% names(cmp$results)))
})

test_that("the baseline pipeline uses the reference model_params and no X", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "selection",
    option_choices = list(list(method = "spearman", top_n = 5)),
    reference_params = list(
      selection_params = list(method = "pearson", top_n = 5),
      model_params      = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  baseline_fit <- cmp$fits[["Baseline"]]
  expect_true(baseline_fit$baseline_model)
  expect_equal(baseline_fit$n_features_input, 0L)
  expect_equal(baseline_fit$model_params$method, "lm")
})

test_that("a failing option is excluded with a message, others still succeed", {
  d <- .make_compare_data()

  msgs <- testthat::capture_messages(
    cmp <- compare_pipelines(
      Y = d$Y, X = d$X,
      option_type    = "selection",
      option_choices = list(
        good = list(method = "spearman", top_n = 5),
        bad  = list(method = "not_a_real_method", top_n = 5)
      ),
      reference_params = list(
        selection_params = list(method = "pearson", top_n = 5),
        model_params      = list(method = "lm")
      ),
      folds   = 4,
      verbose = TRUE
    )
  )

  expect_true(any(grepl("'bad' failed", msgs)))
  expect_false("bad" %in% cmp$results$pipeline)
  expect_true("good" %in% cmp$results$pipeline)
  expect_true("Baseline" %in% cmp$results$pipeline)
  expect_true("Reference" %in% cmp$results$pipeline)
})

test_that("same-method options are auto-labelled with a numeric suffix", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "selection",
    option_choices = list(
      list(method = "spearman", top_n = 5),
      list(method = "spearman", top_n = 10)
    ),
    reference_params = list(
      selection_params = list(method = "pearson", top_n = 5),
      model_params      = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  expect_setequal(names(cmp$option_choices), c("spearman_1", "spearman_2"))
  expect_setequal(cmp$results$pipeline,
                  c("Baseline", "Reference", "spearman_1", "spearman_2"))
})

test_that("user-supplied option names are respected", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "selection",
    option_choices = list(
      loose = list(method = "spearman", top_n = 5),
      tight = list(method = "spearman", top_n = 15)
    ),
    reference_params = list(
      selection_params = list(method = "pearson", top_n = 5),
      model_params      = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  expect_setequal(names(cmp$option_choices), c("loose", "tight"))
})

test_that("gene_level_fc engineering options restrict other pipelines to timepoint == 1", {
  n_ind <- 12
  p     <- 6
  set.seed(2)
  n <- n_ind * 2
  individual_id <- rep(seq_len(n_ind), each = 2)
  timepoint     <- rep(c(0, 1), times = n_ind)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- rnorm(n)

  cmp <- compare_pipelines(
    Y = Y, X = X,
    option_type    = "engineering",
    option_choices = list(
      fc = list(method = "engineer", gene_level_fc = TRUE)
    ),
    reference_params = list(
      engineering_params = list(method = "engineer"),
      model_params        = list(method = "lm")
    ),
    individual_id = individual_id,
    timepoint     = timepoint,
    folds   = 3,
    verbose = FALSE
  )

  expect_equal(length(cmp$fits[["fc"]]$observed), n_ind)
  expect_equal(length(cmp$fits[["Reference"]]$observed), n_ind)
  expect_equal(length(cmp$fits[["Baseline"]]$observed), n_ind)
})

test_that("gene_level_fc engineering option errors informatively without individual_id/timepoint", {
  d <- .make_compare_data()

  expect_error(
    compare_pipelines(
      Y = d$Y, X = d$X,
      option_type    = "engineering",
      option_choices = list(fc = list(method = "engineer", gene_level_fc = TRUE)),
      reference_params = list(
        engineering_params = list(method = "engineer"),
        model_params        = list(method = "lm")
      ),
      folds   = 4,
      verbose = FALSE
    ),
    "individual_id"
  )
})

test_that("individual_id/timepoint are passed through for non-engineering option_type comparisons (regression test)", {
  testthat::skip_if_not_installed("dearseq")

  n_ind <- 10
  p     <- 6
  set.seed(3)
  n <- n_ind * 2
  individual_id <- rep(seq_len(n_ind), each = 2)
  timepoint     <- rep(c(0, 1), times = n_ind)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- rnorm(n)

  # option_type = "selection" (not "engineering") with a dearseq paired-mode
  # option: individual_id/timepoint must still reach predict_cv() for every
  # pipeline, since the gene_level_fc row-parity mechanism (which used to
  # gate this pass-through) is irrelevant here.
  cmp <- compare_pipelines(
    Y = Y, X = X,
    option_type    = "selection",
    option_choices = list(
      dearseq_paired = list(method = "dearseq", dearseq_mode = "paired",
                           threshold = 1.1)
    ),
    reference_params = list(
      selection_params = list(method = "variance", top_n = 3),
      model_params      = list(method = "lm")
    ),
    individual_id = individual_id,
    timepoint     = timepoint,
    folds   = 2,
    verbose = FALSE
  )

  expect_true("dearseq_paired" %in% cmp$results$pipeline)
  expect_true("Baseline" %in% cmp$results$pipeline)
  expect_true("Reference" %in% cmp$results$pipeline)
})

test_that("print and plot methods run without error", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "selection",
    option_choices = list(
      list(method = "spearman", top_n = 5),
      list(method = "variance", top_n = 5)
    ),
    reference_params = list(
      selection_params = list(method = "pearson", top_n = 5),
      model_params      = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  expect_output(print(cmp))
  p1 <- plot(cmp)
  expect_s3_class(p1, "ggplot")
  p2 <- plot(cmp, metric = "all")
  expect_s3_class(p2, "ggplot")
})

test_that("plot.predictomics_comparison shows the baseline as a line, not a bar", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "selection",
    option_choices = list(
      list(method = "spearman", top_n = 5),
      list(method = "variance", top_n = 5)
    ),
    reference_params = list(
      selection_params = list(method = "pearson", top_n = 5),
      model_params      = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  p <- plot(cmp)

  # First layer is the bar geom_col: one bar per Reference/Alternative
  # pipeline, none for the Baseline (which is a dashed line layer instead).
  bar_data <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(bar_data), nrow(cmp$results) - 1L)

  # A linetype scale exists for the "Baseline" reference line, grouped
  # visually under the bar fill scale's "Pipelines" legend title without
  # repeating that title a second time (linetype has no title of its own).
  expect_equal(p$labels$fill, "Pipelines")
  expect_null(p$labels$linetype)

  has_linetype_scale <- any(vapply(
    p$scales$scales, function(s) "linetype" %in% s$aesthetics, logical(1)
  ))
  expect_true(has_linetype_scale)
})

test_that("plot.predictomics_comparison marks the best-performing bar without changing fill colours", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "selection",
    option_choices = list(
      list(method = "spearman", top_n = 5),
      list(method = "variance", top_n = 5)
    ),
    reference_params = list(
      selection_params = list(method = "pearson", top_n = 5),
      model_params      = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  p <- plot(cmp, metric = "sRMSE")

  expect_true("is_best" %in% names(p$data))
  expect_equal(sum(p$data$is_best), 1L)

  best_pipeline     <- as.character(p$data$pipeline[p$data$is_best])
  expected_best     <- as.character(p$data$pipeline[which.min(p$data$metric_value)])
  expect_equal(best_pipeline, expected_best)

  # Only two fills are ever used (Reference/Alternative) - is_best never
  # introduces a third fill level.
  expect_equal(nlevels(droplevels(p$data$role)), 2L)

  # A "pattern" scale exists for the hatched "Best pipeline" highlight,
  # separate from the "Pipelines" fill legend.
  has_pattern_scale <- any(vapply(
    p$scales$scales, function(s) "pattern" %in% s$aesthetics, logical(1)
  ))
  expect_true(has_pattern_scale)
})

test_that("predict_cv's own messages are suppressed, but compare_pipelines's are not", {
  d <- .make_compare_data()

  msgs <- testthat::capture_messages(
    cmp <- compare_pipelines(
      Y = d$Y, X = d$X,
      option_type    = "selection",
      option_choices = list(
        list(method = "spearman", top_n = 5)
      ),
      reference_params = list(
        selection_params = list(method = "pearson", top_n = 5),
        model_params      = list(method = "glmnet")
      ),
      folds   = 4,
      verbose = TRUE
    )
  )

  expect_true(any(grepl("Fitting pipeline", msgs)))
  expect_false(any(grepl("Double variable selection", msgs)))
  expect_false(any(grepl("Starting", msgs)))
})


# -----------------------------------------------------------------------------
# option_type = "predictors" / "response"
# -----------------------------------------------------------------------------

test_that("compare_pipelines runs a basic predictors comparison", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "predictors",
    option_choices = list(
      first_half  = d$X[, 1:10],
      second_half = d$X[, 11:20]
    ),
    reference_params = list(model_params = list(method = "lm")),
    folds   = 4,
    verbose = FALSE
  )

  expect_s3_class(cmp, "predictomics_comparison")
  expect_equal(nrow(cmp$results), 4L)
  expect_setequal(cmp$results$pipeline,
                  c("Baseline", "Reference", "first_half", "second_half"))

  # The reference pipeline used the original (full) X: n_features_input == p
  expect_equal(cmp$fits[["Reference"]]$n_features_input, ncol(d$X))
  expect_equal(cmp$fits[["first_half"]]$n_features_input, 10L)
  expect_equal(cmp$fits[["second_half"]]$n_features_input, 10L)
})

test_that("predictors option respects reference_params$X for the Reference row", {
  d <- .make_compare_data()

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "predictors",
    option_choices = list(alt = d$X[, 1:5]),
    reference_params = list(
      X            = d$X[, 1:8],
      model_params = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  expect_equal(cmp$fits[["Reference"]]$n_features_input, 8L)
  expect_equal(cmp$fits[["alt"]]$n_features_input, 5L)
})

test_that("predictors option errors informatively on wrong nrow", {
  d <- .make_compare_data()

  expect_error(
    compare_pipelines(
      Y = d$Y, X = d$X,
      option_type    = "predictors",
      option_choices = list(bad = d$X[1:5, ]),
      reference_params = list(model_params = list(method = "lm")),
      folds   = 4,
      verbose = FALSE
    ),
    "nrow"
  )
})

test_that("compare_pipelines runs a basic response comparison", {
  d <- .make_compare_data()
  Y_alt <- d$X[, 2] * 3 + rnorm(nrow(d$X))

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "response",
    option_choices = list(alt_response = Y_alt),
    reference_params = list(model_params = list(method = "lm")),
    folds   = 4,
    verbose = FALSE
  )

  expect_s3_class(cmp, "predictomics_comparison")
  expect_setequal(cmp$results$pipeline,
                  c("Baseline", "Reference", "alt_response"))

  expect_equal(cmp$fits[["Reference"]]$observed, d$Y, ignore_attr = TRUE)
  expect_equal(cmp$fits[["alt_response"]]$observed, Y_alt, ignore_attr = TRUE)
  # Baseline is computed against the reference Y, not the alternate response
  expect_equal(cmp$fits[["Baseline"]]$observed, d$Y, ignore_attr = TRUE)
})

test_that("response option respects reference_params$Y for the Reference/Baseline rows", {
  d <- .make_compare_data()
  Y_ref <- d$X[, 3] * 2 + rnorm(nrow(d$X))
  Y_alt <- d$X[, 4] * 2 + rnorm(nrow(d$X))

  cmp <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "response",
    option_choices = list(alt = Y_alt),
    reference_params = list(
      Y            = Y_ref,
      model_params = list(method = "lm")
    ),
    folds   = 4,
    verbose = FALSE
  )

  expect_equal(cmp$fits[["Reference"]]$observed, Y_ref, ignore_attr = TRUE)
  expect_equal(cmp$fits[["Baseline"]]$observed, Y_ref, ignore_attr = TRUE)
  expect_equal(cmp$fits[["alt"]]$observed, Y_alt, ignore_attr = TRUE)
})

test_that("response option errors informatively on wrong length", {
  d <- .make_compare_data()

  expect_error(
    compare_pipelines(
      Y = d$Y, X = d$X,
      option_type    = "response",
      option_choices = list(bad = d$Y[1:5]),
      reference_params = list(model_params = list(method = "lm")),
      folds   = 4,
      verbose = FALSE
    ),
    "length"
  )
})

test_that("unnamed predictors/response options fall back to option_type-based labels", {
  d <- .make_compare_data()
  Y_alt <- d$X[, 2] * 2 + rnorm(nrow(d$X))

  cmp_x <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "predictors",
    option_choices = list(d$X[, 1:5], d$X[, 6:10]),
    reference_params = list(model_params = list(method = "lm")),
    folds   = 4,
    verbose = FALSE
  )
  expect_setequal(names(cmp_x$option_choices), c("predictors_1", "predictors_2"))

  cmp_y <- compare_pipelines(
    Y = d$Y, X = d$X,
    option_type    = "response",
    option_choices = list(Y_alt),
    reference_params = list(model_params = list(method = "lm")),
    folds   = 4,
    verbose = FALSE
  )
  expect_setequal(names(cmp_y$option_choices), "response")
})
