test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

# -----------------------------------------------------------------------------
# gene_level_fc engineering
# -----------------------------------------------------------------------------

.make_gene_level_fc_data <- function(n_ind = 10, p = 5, seed = 1) {
  set.seed(seed)
  n <- n_ind * 2
  individual_id <- rep(seq_len(n_ind), each = 2)
  timepoint     <- rep(c(0, 1), times = n_ind)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- rnorm(n)
  list(X = X, Y = Y, individual_id = individual_id, timepoint = timepoint)
}

test_that("gene_level_fc collapses paired rows into post - pre fold changes", {
  d <- .make_gene_level_fc_data()

  result <- predict_cv(
    Y = d$Y, X = d$X,
    folds = 2, cv_type = "kfold",
    engineering_params = list(method = "engineer", gene_level_fc = TRUE),
    individual_id = d$individual_id,
    timepoint = d$timepoint,
    verbose = FALSE
  )

  expect_equal(length(result$observed), 10)
  expect_equal(result$n_samples, 10)
})

test_that(".compute_gene_level_fc computes post - pre correctly", {
  d <- .make_gene_level_fc_data(n_ind = 4, p = 3)
  out <- .compute_gene_level_fc(d$X, d$Y, d$individual_id, d$timepoint)

  for (i in seq_len(4)) {
    idx <- which(d$individual_id == i)
    pre  <- idx[d$timepoint[idx] == 0]
    post <- idx[d$timepoint[idx] == 1]
    expect_equal(out$X[i, ], d$X[post, ] - d$X[pre, ], ignore_attr = TRUE)
    expect_equal(out$Y[i], d$Y[post])
  }
})

test_that("gene_level_fc errors when individual_id/timepoint missing", {
  d <- .make_gene_level_fc_data()

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", gene_level_fc = TRUE),
      verbose = FALSE
    ),
    "individual_id"
  )
})

test_that("gene_level_fc errors when an individual lacks exactly two observations", {
  d <- .make_gene_level_fc_data()
  individual_id <- d$individual_id
  individual_id[1] <- 999  # orphan a single pre-treatment row

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", gene_level_fc = TRUE),
      individual_id = individual_id,
      timepoint = d$timepoint,
      verbose = FALSE
    ),
    "exactly 2 observations"
  )
})

test_that("gene_level_fc errors when an individual has two obs at the same timepoint", {
  d <- .make_gene_level_fc_data()
  timepoint <- d$timepoint
  timepoint[2] <- 0  # individual 1 now has two pre-treatment rows

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", gene_level_fc = TRUE),
      individual_id = d$individual_id,
      timepoint = timepoint,
      verbose = FALSE
    ),
    "one pre-treatment"
  )
})

test_that("gene_level_fc is incompatible with rise_paired = TRUE", {
  d <- .make_gene_level_fc_data()

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      engineering_params = list(method = "engineer", gene_level_fc = TRUE),
      selection_params   = list(method = "rise", rise_paired = TRUE,
                                rise_epsilon = 0.1),
      individual_id = d$individual_id,
      timepoint = d$timepoint,
      treatment = rep(c(0, 1), 10),
      verbose = FALSE
    ),
    "not compatible"
  )
})

test_that("gene_level_fc allowed when rise_paired = FALSE", {
  d <- .make_gene_level_fc_data()

  expect_error(
    predict_cv(
      Y = d$Y, X = d$X,
      folds = 2,
      engineering_params = list(method = "engineer", gene_level_fc = TRUE),
      selection_params   = list(method = "variance", top_n = 3),
      individual_id = d$individual_id,
      timepoint = d$timepoint,
      verbose = FALSE
    ),
    NA
  )
})
