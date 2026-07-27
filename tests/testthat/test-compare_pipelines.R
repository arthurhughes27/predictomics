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
