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
  genesets <- list(
    setA = paste0("gene", 1:4),
    setB = paste0("gene", 5:8)
  )
  list(X = X, Y = Y, individual_id = individual_id, timepoint = timepoint,
       treatment = treatment, genesets = genesets)
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

test_that(".build_dearseq_covariates drops a single-level covariate with a warning", {
  covs <- data.frame(sex = rep("Female", 6), age = rnorm(6))

  expect_warning(design <- .build_dearseq_covariates(covs), "sex")
  expect_false("sexFemale" %in% colnames(design))
  expect_true("age" %in% colnames(design))
})

test_that(".build_dearseq_covariates returns NULL with a message when all covariates are single-level", {
  covs <- data.frame(sex = rep("Female", 6))

  expect_message(design <- .build_dearseq_covariates(covs),
                 "All dearseq covariates")
  expect_null(design)
})

test_that("run_selection no longer errors on a constant covariate for dearseq (regression test)", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()
  constant_covariates <- data.frame(sex = rep("Female", nrow(d$X)))

  expect_warning(
    res <- run_selection(
      X_train    = d$X,
      covariates = constant_covariates,
      treatment  = d$treatment,
      params     = list(method = "dearseq", dearseq_mode = "classic",
                        threshold = 1.1)
    ),
    "sex"
  )
  expect_equal(length(res$selection_scores), ncol(d$X))
})

test_that("run_selection validates selection_params$dearseq_covariates", {
  d <- .make_dearseq_data()

  expect_error(
    run_selection(
      X_train    = d$X,
      treatment  = d$treatment,
      params     = list(method = "dearseq", dearseq_mode = "classic",
                        threshold = 0.05,
                        dearseq_covariates = matrix(1:3, nrow = 3))
    ),
    "same number of rows"
  )

  bad_names <- matrix(rnorm(nrow(d$X) * 2), nrow = nrow(d$X))
  expect_error(
    run_selection(
      X_train    = d$X,
      treatment  = d$treatment,
      params     = list(method = "dearseq", dearseq_mode = "classic",
                        threshold = 0.05,
                        dearseq_covariates = bad_names)
    ),
    "column names"
  )

  expect_error(
    run_selection(
      X_train    = d$X,
      treatment  = d$treatment,
      params     = list(method = "dearseq", dearseq_mode = "classic",
                        threshold = 0.05,
                        dearseq_covariates = 1:nrow(d$X))
    ),
    "matrix or data frame"
  )
})

test_that("run_selection uses dearseq_covariates in place of covariates for the dearseq comparison", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  # A covariate constant across every row would otherwise be dropped (with a
  # warning) by .build_dearseq_covariates() - supplying a non-constant
  # dearseq_covariates instead should run cleanly, with no warning about it.
  constant_covariates <- data.frame(sex = rep("Female", nrow(d$X)))
  alt_covariates       <- data.frame(age = rnorm(nrow(d$X)))

  res <- run_selection(
    X_train            = d$X,
    covariates         = constant_covariates,
    treatment          = d$treatment,
    params             = list(method = "dearseq", dearseq_mode = "classic",
                              threshold = 1.1,
                              dearseq_covariates = alt_covariates)
  )

  expect_equal(length(res$selection_scores), ncol(d$X))
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

# -----------------------------------------------------------------------------
# dearseq_level = "geneset"
# -----------------------------------------------------------------------------

test_that("dearseq_level defaults to 'gene' and 'geneset' requires genesets", {
  d <- .make_dearseq_data()

  expect_error(
    .validate_selection_params(
      list(method = "dearseq", threshold = 0.05, dearseq_level = "geneset"),
      p = ncol(d$X)
    ),
    "genesets"
  )

  expect_error(
    .validate_selection_params(
      list(method = "dearseq", threshold = 0.05, dearseq_level = "geneset",
           genesets = d$genesets),
      p = ncol(d$X)
    ),
    NA
  )
})

test_that("dearseq_level must be 'gene' or 'geneset'", {
  d <- .make_dearseq_data()

  expect_error(
    .validate_selection_params(
      list(method = "dearseq", threshold = 0.05, dearseq_level = "bad"),
      p = ncol(d$X)
    ),
    "dearseq_level"
  )
})

test_that("top_n for dearseq geneset-level is not bounded at validation time (checked by run_selection instead)", {
  d <- .make_dearseq_data()

  # top_n (3) exceeds the number of genesets (2), but this is no longer a
  # validation error - run_selection() handles it directly (see below).
  expect_error(
    .validate_selection_params(
      list(method = "dearseq", dearseq_level = "geneset", genesets = d$genesets,
           top_n = 3),
      p = ncol(d$X)
    ),
    NA
  )

  expect_error(
    .validate_selection_params(
      list(method = "dearseq", dearseq_level = "geneset", genesets = d$genesets,
           top_n = 2),
      p = ncol(d$X)
    ),
    NA
  )
})

test_that("run_selection selects all genesets (with a message) when top_n exceeds the number of genesets", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  msgs <- testthat::capture_messages(
    res <- run_selection(
      X_train       = d$X,
      individual_id = d$individual_id,
      timepoint     = d$timepoint,
      params        = list(method = "dearseq", dearseq_mode = "paired",
                           dearseq_level = "geneset", genesets = d$genesets,
                           top_n = 3)
    )
  )

  expect_true(any(grepl("exceeds \\(or equals\\)", msgs)))
  expect_setequal(res$selected_features, unlist(d$genesets, use.names = FALSE))
  expect_true(all(is.na(res$selection_scores)))
})

test_that("run_selection resolves selected genesets to member gene names (geneset-level)", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  res <- run_selection(
    X_train   = d$X,
    treatment = d$treatment,
    params    = list(method = "dearseq", dearseq_mode = "classic",
                     dearseq_level = "geneset", genesets = d$genesets,
                     threshold = 1.1)
  )

  expect_equal(length(res$selection_scores), length(d$genesets))
  expect_setequal(names(res$selection_scores), names(d$genesets))
  expect_true(all(res$selected_features %in% colnames(d$X)))
  expect_setequal(res$selected_features, colnames(d$X))
})

test_that("predict_cv stops when dearseq_level = 'gene' is combined with engineering genesets", {
  d <- .make_dearseq_data()

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      selection_params   = list(method = "dearseq", dearseq_mode = "classic",
                                threshold = 0.05),
      engineering_params = list(method = "engineer", genesets = d$genesets,
                                agg_method = "mean"),
      treatment = d$treatment,
      verbose   = FALSE
    ),
    "not compatible"
  )
})

test_that("predict_cv allows dearseq_level = 'geneset' combined with engineering genesets", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      folds = 2,
      selection_params   = list(method = "dearseq", dearseq_mode = "classic",
                                dearseq_level = "geneset", genesets = d$genesets,
                                threshold = 1.1),
      engineering_params = list(method = "engineer", genesets = d$genesets,
                                agg_method = "mean"),
      treatment = d$treatment,
      verbose   = FALSE
    ),
    NA
  )
})

test_that("predict_cv drops engineering genesets with no surviving genes after geneset-level dearseq filtering", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  msgs <- testthat::capture_messages(
    result <- predict_cv(
      Y = d$Y, X = d$X,
      folds = 2,
      selection_params   = list(method = "dearseq", dearseq_mode = "classic",
                                dearseq_level = "geneset", genesets = d$genesets,
                                top_n = 1),
      engineering_params = list(method = "engineer", genesets = d$genesets,
                                agg_method = "mean"),
      treatment = d$treatment,
      verbose   = TRUE
    )
  )

  # top_n = 1 with two disjoint, equal-sized genesets guarantees exactly one
  # geneset survives, leaving the other with zero overlap in X; predict_cv
  # must drop it automatically rather than erroring downstream.
  expect_true(any(grepl("Removed 1 geneset", msgs)))
  expect_equal(result$dearseq_selection$n_selected, 1)
  expect_true(result$dearseq_selection$selected_features %in% names(d$genesets))
  expect_false(is.null(result$predicted))
})

test_that("dearseq_selection reports geneset names (not gene names) when dearseq_level = 'geneset'", {
  testthat::skip_if_not_installed("dearseq")
  d <- .make_dearseq_data()

  result <- predict_cv(
    Y = d$Y, X = d$X,
    folds = 2,
    selection_params = list(method = "dearseq", dearseq_mode = "classic",
                            dearseq_level = "geneset", genesets = d$genesets,
                            top_n = 1),
    treatment = d$treatment,
    verbose   = FALSE
  )

  expect_equal(result$dearseq_selection$n_selected, 1)
  expect_length(result$dearseq_selection$selected_features, 1)
  expect_true(result$dearseq_selection$selected_features %in% names(d$genesets))
  expect_setequal(names(result$dearseq_selection$selection_scores), names(d$genesets))
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
