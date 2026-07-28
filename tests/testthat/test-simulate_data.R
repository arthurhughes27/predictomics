# -----------------------------------------------------------------------------
# simulate_predictomics_data()
# -----------------------------------------------------------------------------

test_that("cross_sectional design returns correctly shaped, complete data", {
  d <- simulate_predictomics_data(n = 30, p = 60, n_genesets = 6,
                                  geneset_size = 10, seed = 1)

  expect_equal(dim(d$X), c(30, 60))
  expect_equal(length(d$Y), 30)
  expect_false(anyNA(d$X))
  expect_false(anyNA(d$Y))
  expect_equal(length(d$treatment), 30)
  expect_true(all(d$treatment %in% c(0, 1)))
  expect_equal(nrow(d$covariates), 30)
  expect_setequal(colnames(d$covariates), c("age", "sex"))
  expect_equal(length(d$genesets), 6)
  expect_equal(length(d$genesets[[1]]), 10)
  expect_setequal(unlist(d$genesets, use.names = FALSE), colnames(d$X))
  expect_true(all(d$signal_genesets %in% names(d$genesets)))
  expect_null(d$individual_id)
})

test_that("paired design returns correctly shaped, complete data with valid pairing", {
  d <- simulate_predictomics_data(n = 20, p = 40, n_genesets = 4,
                                  geneset_size = 10, design = "paired", seed = 1)

  expect_equal(dim(d$X), c(40, 40))
  expect_equal(length(d$Y), 40)
  expect_false(anyNA(d$X))
  expect_false(anyNA(d$Y))
  expect_equal(length(d$individual_id), 40)
  expect_equal(length(d$timepoint), 40)
  expect_true(all(table(d$individual_id) == 2))

  for (id in unique(d$individual_id)) {
    tp <- d$timepoint[d$individual_id == id]
    expect_setequal(tp, c(0, 1))
  }

  expect_equal(nrow(d$covariates), 40)
})

test_that("simulated signal genes are more predictive than noise genes", {
  d <- simulate_predictomics_data(n = 100, p = 100, n_genesets = 10,
                                  geneset_size = 10, seed = 1)

  signal_genes <- unlist(d$genesets[d$signal_genesets])
  noise_genes  <- setdiff(colnames(d$X), signal_genes)

  signal_cor <- mean(abs(apply(d$X[, signal_genes], 2, cor, y = d$Y)))
  noise_cor  <- mean(abs(apply(d$X[, noise_genes], 2, cor, y = d$Y)))

  expect_gt(signal_cor, noise_cor)
})

test_that("design must be a supported value", {
  expect_error(
    simulate_predictomics_data(design = "not_a_design"),
    "design"
  )
})

test_that("p must equal n_genesets * geneset_size", {
  expect_error(
    simulate_predictomics_data(p = 99, n_genesets = 10, geneset_size = 20),
    "n_genesets \\* geneset_size"
  )
})

test_that("predict_cv runs end-to-end on simulated cross-sectional data", {
  d <- simulate_predictomics_data(n = 30, p = 40, n_genesets = 4,
                                  geneset_size = 10, seed = 1)

  result <- predict_cv(
    Y = d$Y, X = d$X,
    selection_params = list(method = "pearson", top_n = 10),
    folds   = 3,
    verbose = FALSE
  )

  expect_s3_class(result, "predictomics")
  expect_equal(length(result$observed), 30)
})

test_that("predict_cv runs end-to-end on simulated paired data with gene_level_fc", {
  d <- simulate_predictomics_data(n = 16, p = 40, n_genesets = 4,
                                  geneset_size = 10, design = "paired", seed = 1)

  result <- predict_cv(
    Y = d$Y, X = d$X,
    engineering_params = list(method = "engineer", gene_level_fc = TRUE),
    selection_params   = list(method = "pearson", top_n = 10),
    individual_id      = d$individual_id,
    timepoint          = d$timepoint,
    folds   = 4,
    verbose = FALSE
  )

  expect_s3_class(result, "predictomics")
  expect_equal(length(result$observed), 16)
})
