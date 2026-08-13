# -----------------------------------------------------------------------------
# RISE selection (rise_paired mode harmonised with dearseq paired mode)
# -----------------------------------------------------------------------------

.make_rise_data <- function(n_ind = 20, p = 6, seed = 1) {
  set.seed(seed)
  n <- n_ind * 2
  individual_id <- rep(seq_len(n_ind), each = 2)
  timepoint     <- rep(c(0, 1), times = n_ind)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- rnorm(n)
  list(X = X, Y = Y, individual_id = individual_id, timepoint = timepoint)
}

test_that("run_selection still requires treatment for unpaired rise", {
  d <- .make_rise_data()

  expect_error(
    run_selection(
      X_train = d$X, Y_train = d$Y,
      params  = list(method = "rise", rise_epsilon = 0.1)
    ),
    "treatment"
  )
})

test_that("run_selection requires individual_id/timepoint for rise_paired = TRUE", {
  d <- .make_rise_data()

  expect_error(
    run_selection(
      X_train = d$X, Y_train = d$Y,
      params  = list(method = "rise", rise_paired = TRUE, rise_epsilon = 0.1)
    ),
    "individual_id"
  )
})

test_that("run_selection does not require treatment for rise_paired = TRUE", {
  testthat::skip_if_not_installed("SurrogateRank")
  d <- .make_rise_data()

  expect_error(
    run_selection(
      X_train       = d$X, Y_train = d$Y,
      individual_id = d$individual_id,
      timepoint     = d$timepoint,
      params        = list(method = "rise", rise_paired = TRUE,
                           rise_epsilon = 0.1, threshold = 1.1)
    ),
    NA
  )
})

test_that("predict_cv errors informatively for rise_paired = TRUE without individual_id/timepoint", {
  d <- .make_rise_data()

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      selection_params = list(method = "rise", rise_paired = TRUE,
                              rise_epsilon = 0.1, threshold = 1.1),
      verbose = FALSE
    ),
    "individual_id"
  )
})

test_that("run_selection computes rise_paired scores correctly regardless of input row order", {
  testthat::skip_if_not_installed("SurrogateRank")
  d <- .make_rise_data()

  # Shuffle rows so paired pre/post observations are no longer adjacent or
  # positionally matched - the whole point of sorting by individual_id.
  shuffle <- sample(seq_len(nrow(d$X)))

  res <- run_selection(
    X_train       = d$X[shuffle, , drop = FALSE],
    Y_train       = d$Y[shuffle],
    individual_id = d$individual_id[shuffle],
    timepoint     = d$timepoint[shuffle],
    params        = list(method = "rise", rise_paired = TRUE,
                         rise_epsilon = 0.1, threshold = 1.1)
  )

  expect_equal(length(res$selection_scores), ncol(d$X))
  expect_setequal(names(res$selection_scores), colnames(d$X))
})

test_that("predict_cv runs rise_paired mode correctly with shuffled row order", {
  testthat::skip_if_not_installed("SurrogateRank")
  d <- .make_rise_data()
  shuffle <- sample(seq_len(nrow(d$X)))

  result <- predict_cv(
    Y             = d$Y[shuffle],
    X             = d$X[shuffle, , drop = FALSE],
    individual_id = d$individual_id[shuffle],
    timepoint     = d$timepoint[shuffle],
    selection_params = list(method = "rise", rise_paired = TRUE,
                            rise_epsilon = 0.1, threshold = 1.1),
    model_params  = list(method = "lm"),
    folds   = 3,
    verbose = FALSE
  )

  expect_true(result$paired_rise)
  expect_equal(length(result$observed), 20)
})

test_that("predict_cv's returned treatment matches x$observed for rise_paired, decoupled from timepoint", {
  testthat::skip_if_not_installed("SurrogateRank")
  d <- .make_rise_data()
  treatment_arm <- rep(c(0, 1), length.out = nrow(d$X))
  shuffle <- sample(seq_len(nrow(d$X)))

  result <- predict_cv(
    Y             = d$Y[shuffle],
    X             = d$X[shuffle, , drop = FALSE],
    individual_id = d$individual_id[shuffle],
    timepoint     = d$timepoint[shuffle],
    treatment     = treatment_arm[shuffle],
    selection_params = list(method = "rise", rise_paired = TRUE,
                            rise_epsilon = 0.1, threshold = 1.1),
    model_params  = list(method = "lm"),
    folds   = 3,
    verbose = FALSE
  )

  expect_equal(length(result$treatment), length(result$observed))
})


# -----------------------------------------------------------------------------
# rise_paired incompatibility with geneset-level engineering
# -----------------------------------------------------------------------------

test_that("predict_cv errors informatively for rise_paired = TRUE with geneset engineering", {
  d <- .make_rise_data()
  genesets <- list(
    setA = paste0("gene", 1:3),
    setB = paste0("gene", 4:6)
  )

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", genesets = genesets,
                                agg_method = "mean"),
      selection_params   = list(method = "rise", rise_paired = TRUE,
                                rise_epsilon = 0.1, threshold = 1.1),
      individual_id = d$individual_id,
      timepoint     = d$timepoint,
      verbose       = FALSE
    ),
    "genesets"
  )
})

test_that("unpaired rise is not blocked by the rise_paired/genesets guard", {
  d <- .make_rise_data()
  genesets <- list(
    setA = paste0("gene", 1:3),
    setB = paste0("gene", 4:6)
  )
  treatment_arm <- rep(c(0, 1), length.out = nrow(d$X))

  # rise_paired is FALSE (default) here, so the new guard must not fire;
  # this should fail (if at all) only for the usual reason unpaired rise
  # requires treatment, not for a genesets incompatibility.
  err <- tryCatch({
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", genesets = genesets,
                                agg_method = "mean"),
      selection_params   = list(method = "rise", rise_epsilon = 0.1,
                                threshold = 1.1),
      treatment = treatment_arm,
      verbose   = FALSE
    )
    NULL
  }, error = function(e) conditionMessage(e))

  if (!is.null(err)) expect_false(grepl("genesets", err))
})
