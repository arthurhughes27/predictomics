# -----------------------------------------------------------------------------
# selection_params$method = "variance" combined with
# engineering_params$col_transform = "z": scores are computed on a separately
# re-engineered, pre-z-score matrix (same aggregation, if any), while
# modelling still uses the fully engineered (z-scored) matrix.
# -----------------------------------------------------------------------------

.make_variance_zscore_data <- function(n = 40, seed = 1) {
  set.seed(seed)
  # Genes 1-3 have much larger raw variance than genes 4-10; z-scoring erases
  # this gap entirely (every column becomes variance ~1), so post-z variance
  # selection could not reliably recover it, but pre-z variance selection
  # should pick genes 1-3 in essentially every fold.
  X <- cbind(
    matrix(rnorm(n * 3, sd = 10), nrow = n, ncol = 3),
    matrix(rnorm(n * 7, sd = 0.1), nrow = n, ncol = 7)
  )
  colnames(X) <- paste0("gene", seq_len(10))
  Y <- X[, 1] * 0.1 + rnorm(n)
  list(X = X, Y = Y)
}

test_that("variance + col_transform = 'z' warns", {
  d <- .make_variance_zscore_data()

  expect_warning(
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", col_transform = "z"),
      selection_params   = list(method = "variance", top_n = 3),
      folds   = 4,
      verbose = FALSE
    ),
    "variance"
  )
})

test_that("variance + col_transform = 'z' selects by pre-z-score variance (outside_cv, no aggregation)", {
  d <- .make_variance_zscore_data()

  result <- suppressWarnings(predict_cv(
    Y = d$Y, X = d$X,
    engineering_params = list(method = "engineer", col_transform = "z"),
    selection_params   = list(method = "variance", top_n = 3),
    outside_cv = TRUE,
    folds   = 4,
    verbose = FALSE
  ))

  expect_setequal(result$outside_cv_selection$selected_features,
                  c("gene1", "gene2", "gene3"))
})

test_that("variance + col_transform = 'z' still models on the z-scored matrix (predictions match an equivalent pre-filtered pipeline)", {
  d <- .make_variance_zscore_data()

  result_auto <- suppressWarnings(predict_cv(
    Y = d$Y, X = d$X,
    engineering_params = list(method = "engineer", col_transform = "z"),
    selection_params   = list(method = "variance", top_n = 3),
    model_params = list(method = "lm"),
    folds   = 4,
    seed    = 999,
    verbose = FALSE
  ))

  # Every fold should have selected exactly genes 1-3, given the large
  # variance gap constructed above.
  for (k in seq_along(result_auto$fold_selection_diagnostics)) {
    expect_setequal(result_auto$fold_selection_diagnostics[[k]]$selected_features,
                    c("gene1", "gene2", "gene3"))
  }

  # An equivalent pipeline that pre-filters to those same 3 genes and applies
  # col_transform = "z" with no selection step should produce identical
  # predictions, confirming modelling used the fully z-scored matrix.
  result_manual <- predict_cv(
    Y = d$Y, X = d$X[, c("gene1", "gene2", "gene3")],
    engineering_params = list(method = "engineer", col_transform = "z"),
    selection_params   = NULL,
    model_params = list(method = "lm"),
    folds   = 4,
    seed    = 999,
    verbose = FALSE
  )

  expect_equal(result_auto$predicted, result_manual$predicted)
})

test_that("variance + col_transform = 'z' with geneset aggregation selects by pre-z-score aggregated variance", {
  set.seed(2)
  n <- 40
  # setA's members are highly correlated (a coherent, high-post-z-variance
  # module even after standardisation) but have small raw variance; setB's
  # members are large-raw-variance but mutually independent. Pre-z-score
  # aggregated variance should favour setB (raw scale dominates); this is
  # the behaviour being tested, not a claim about which is "right".
  base_signal <- rnorm(n)
  setA <- sapply(1:4, function(i) 0.1 * base_signal + rnorm(n, sd = 0.01))
  setB <- matrix(rnorm(n * 4, sd = 10), nrow = n, ncol = 4)
  X <- cbind(setA, setB)
  colnames(X) <- c(paste0("geneA", 1:4), paste0("geneB", 1:4))
  Y <- rnorm(n)

  genesets <- list(setA = paste0("geneA", 1:4), setB = paste0("geneB", 1:4))

  result <- suppressWarnings(predict_cv(
    Y = Y, X = X,
    engineering_params = list(method = "engineer", col_transform = "z",
                              genesets = genesets, agg_method = "mean"),
    selection_params   = list(method = "variance", top_n = 1),
    outside_cv = TRUE,
    folds   = 4,
    verbose = FALSE
  ))

  expect_equal(result$outside_cv_selection$selected_features, "setB")

  # Sanity check: the same setup without col_transform = "z" (aggregation
  # only) selects on the identical (already pre-z-score) aggregated matrix,
  # so it should agree with the z-score case above.
  result_no_z <- suppressWarnings(predict_cv(
    Y = Y, X = X,
    engineering_params = list(method = "engineer", genesets = genesets,
                              agg_method = "mean"),
    selection_params   = list(method = "variance", top_n = 1),
    outside_cv = TRUE,
    folds   = 4,
    verbose = FALSE
  ))
  expect_equal(result$outside_cv_selection$selected_features,
              result_no_z$outside_cv_selection$selected_features)
})

test_that("pearson/spearman selection is unaffected by col_transform = 'z' (no warning)", {
  d <- .make_variance_zscore_data()

  for (m in c("pearson", "spearman")) {
    expect_warning(
      predict_cv(
        Y = d$Y, X = d$X,
        engineering_params = list(method = "engineer", col_transform = "z"),
        selection_params   = list(method = m, top_n = 3),
        folds   = 4,
        verbose = FALSE
      ),
      NA
    )
  }
})

test_that("variance selection without col_transform = 'z' is unaffected (no warning)", {
  d <- .make_variance_zscore_data()

  expect_warning(
    predict_cv(
      Y = d$Y, X = d$X,
      selection_params = list(method = "variance", top_n = 3),
      folds   = 4,
      verbose = FALSE
    ),
    NA
  )

  expect_warning(
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", col_transform = "none"),
      selection_params   = list(method = "variance", top_n = 3),
      folds   = 4,
      verbose = FALSE
    ),
    NA
  )
})
