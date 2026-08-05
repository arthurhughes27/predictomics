# -----------------------------------------------------------------------------
# ssgsea aggregation
# -----------------------------------------------------------------------------

.make_ssgsea_data <- function(n = 20, p = 10, seed = 1) {
  set.seed(seed)
  X <- matrix(rexp(n * p, rate = 1), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  rownames(X) <- paste0("sample", seq_len(n))
  genesets <- list(
    setA = paste0("gene", 1:5),
    setB = paste0("gene", 6:10)
  )
  list(X = X, genesets = genesets)
}

test_that("engineering_params$agg_method = 'ssgsea' is a valid option", {
  d <- .make_ssgsea_data()
  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "ssgsea")
  expect_error(.validate_engineering_params(params), NA)
})

test_that("ssgsea-specific params are validated", {
  d <- .make_ssgsea_data()

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "ssgsea",
      ssgsea_alpha = -1
    )),
    "ssgsea_alpha"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "ssgsea",
      ssgsea_min_size = 0
    )),
    "ssgsea_min_size"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "ssgsea",
      ssgsea_max_size = 1, ssgsea_min_size = 5
    )),
    "ssgsea_max_size"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "ssgsea",
      ssgsea_normalize = "yes"
    )),
    "ssgsea_normalize"
  )
})

test_that("run_engineering computes ssgsea scores with correct shape and names", {
  testthat::skip_if_not_installed("GSVA")
  d <- .make_ssgsea_data()

  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "ssgsea")
  fit <- run_engineering(X_train = d$X, params = params)

  expect_equal(dim(fit$X_transformed), c(nrow(d$X), length(d$genesets)))
  expect_equal(colnames(fit$X_transformed), names(d$genesets))
  expect_equal(rownames(fit$X_transformed), rownames(d$X))
})

test_that("predict_engineering recomputes ssgsea scores on new data", {
  testthat::skip_if_not_installed("GSVA")
  d <- .make_ssgsea_data(n = 20, seed = 1)
  d_test <- .make_ssgsea_data(n = 8, seed = 2)

  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "ssgsea")
  fit <- run_engineering(X_train = d$X, params = params)

  X_test_transformed <- predict_engineering(fit$fit, X_new = d_test$X)

  expect_equal(dim(X_test_transformed), c(nrow(d_test$X), length(d$genesets)))
  expect_equal(colnames(X_test_transformed), names(d$genesets))
})

test_that("ssgsea errors informatively when a geneset has no overlapping features", {
  testthat::skip_if_not_installed("GSVA")
  d <- .make_ssgsea_data()
  genesets <- d$genesets
  genesets$setC <- c("not_a_gene1", "not_a_gene2")

  params <- list(method = "engineer", col_transform = "none",
                 genesets = genesets, agg_method = "ssgsea")

  expect_error(
    run_engineering(X_train = d$X, params = params),
    "setC"
  )
})

test_that("using col_transform = 'z' with agg_method = 'ssgsea' warns", {
  testthat::skip_if_not_installed("GSVA")
  d <- .make_ssgsea_data()
  params <- list(method = "engineer", col_transform = "z",
                 genesets = d$genesets, agg_method = "ssgsea")

  expect_warning(
    run_engineering(X_train = d$X, params = params),
    "col_transform"
  )
})


# -----------------------------------------------------------------------------
# gsva aggregation
# -----------------------------------------------------------------------------

test_that("engineering_params$agg_method = 'gsva' is a valid option", {
  d <- .make_ssgsea_data()
  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "gsva")
  expect_error(.validate_engineering_params(params), NA)
})

test_that("gsva-specific params are validated", {
  d <- .make_ssgsea_data()

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "gsva",
      gsva_kcdf = "not_a_kernel"
    )),
    "gsva_kcdf"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "gsva",
      gsva_tau = -1
    )),
    "gsva_tau"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "gsva",
      gsva_max_diff = "yes"
    )),
    "gsva_max_diff"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "gsva",
      gsva_abs_ranking = "yes"
    )),
    "gsva_abs_ranking"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "gsva",
      gsva_min_size = 0
    )),
    "gsva_min_size"
  )

  expect_error(
    .validate_engineering_params(list(
      method = "engineer", genesets = d$genesets, agg_method = "gsva",
      gsva_max_size = 1, gsva_min_size = 5
    )),
    "gsva_max_size"
  )
})

test_that("run_engineering computes gsva scores with correct shape and names", {
  testthat::skip_if_not_installed("GSVA")
  d <- .make_ssgsea_data()

  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "gsva")
  fit <- run_engineering(X_train = d$X, params = params)

  expect_equal(dim(fit$X_transformed), c(nrow(d$X), length(d$genesets)))
  expect_equal(colnames(fit$X_transformed), names(d$genesets))
  expect_equal(rownames(fit$X_transformed), rownames(d$X))
})

test_that("predict_engineering recomputes gsva scores on new data", {
  testthat::skip_if_not_installed("GSVA")
  d <- .make_ssgsea_data(n = 20, seed = 1)
  d_test <- .make_ssgsea_data(n = 8, seed = 2)

  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "gsva")
  fit <- run_engineering(X_train = d$X, params = params)

  X_test_transformed <- predict_engineering(fit$fit, X_new = d_test$X)

  expect_equal(dim(X_test_transformed), c(nrow(d_test$X), length(d$genesets)))
  expect_equal(colnames(X_test_transformed), names(d$genesets))
})

test_that("gsva errors informatively when a geneset has no overlapping features", {
  testthat::skip_if_not_installed("GSVA")
  d <- .make_ssgsea_data()
  genesets <- d$genesets
  genesets$setC <- c("not_a_gene1", "not_a_gene2")

  params <- list(method = "engineer", col_transform = "none",
                 genesets = genesets, agg_method = "gsva")

  expect_error(
    run_engineering(X_train = d$X, params = params),
    "setC"
  )
})

# -----------------------------------------------------------------------------
# max aggregation
# -----------------------------------------------------------------------------

test_that("engineering_params$agg_method = 'max' is a valid option", {
  d <- .make_ssgsea_data()
  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "max")
  expect_error(.validate_engineering_params(params), NA)
})

test_that("run_engineering computes the row-wise max within each geneset", {
  d <- .make_ssgsea_data()

  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "max")
  fit <- run_engineering(X_train = d$X, params = params)

  expect_equal(dim(fit$X_transformed), c(nrow(d$X), length(d$genesets)))
  expect_equal(colnames(fit$X_transformed), names(d$genesets))
  expect_equal(fit$X_transformed[, "setA"],
              apply(d$X[, d$genesets$setA], 1, max), ignore_attr = TRUE)
  expect_equal(fit$X_transformed[, "setB"],
              apply(d$X[, d$genesets$setB], 1, max), ignore_attr = TRUE)
})

test_that("predict_engineering recomputes the max on new data", {
  d      <- .make_ssgsea_data(n = 20, seed = 1)
  d_test <- .make_ssgsea_data(n = 8, seed = 2)

  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "max")
  fit <- run_engineering(X_train = d$X, params = params)

  X_test_transformed <- predict_engineering(fit$fit, X_new = d_test$X)

  expect_equal(dim(X_test_transformed), c(nrow(d_test$X), length(d$genesets)))
  expect_equal(X_test_transformed[, "setA"],
              apply(d_test$X[, d$genesets$setA], 1, max), ignore_attr = TRUE)
})


# -----------------------------------------------------------------------------
# NA support
# -----------------------------------------------------------------------------

test_that("col_transform = 'z' tolerates NA via na.rm", {
  d <- .make_ssgsea_data()
  X_na <- d$X
  X_na[1, 1] <- NA
  X_na[2, 2] <- NA

  params <- list(method = "engineer", col_transform = "z")
  fit <- run_engineering(X_train = X_na, params = params)

  expect_false(anyNA(fit$fit$col_means))
  expect_false(anyNA(fit$fit$col_sds))
  expect_true(anyNA(fit$X_transformed))  # the original NA entries remain NA
})

test_that("col_transform = 'z' warns and skips a degenerate (< 2 non-missing) column", {
  d <- .make_ssgsea_data(n = 10)
  X_na <- d$X
  X_na[2:10, 1] <- NA  # only 1 non-missing value in column 1

  params <- list(method = "engineer", col_transform = "z")
  expect_warning(
    fit <- run_engineering(X_train = X_na, params = params),
    "fewer than 2 non-missing"
  )
  expect_equal(fit$fit$col_sds[["gene1"]], 1)
})

test_that("mean/median/sum aggregation tolerate NA via na.rm", {
  d <- .make_ssgsea_data()
  X_na <- d$X
  X_na[1, "gene1"] <- NA

  for (m in c("mean", "median", "sum")) {
    params <- list(method = "engineer", col_transform = "none",
                   genesets = d$genesets, agg_method = m)
    fit <- run_engineering(X_train = X_na, params = params)
    expect_false(anyNA(fit$X_transformed[-1, "setA"]))
    expect_false(is.na(fit$X_transformed[1, "setA"]))  # only 1/5 genes NA
  }
})

test_that("aggregation yields NA when all geneset members are NA for a sample", {
  d <- .make_ssgsea_data()
  X_na <- d$X
  X_na[1, d$genesets$setA] <- NA

  params <- list(method = "engineer", col_transform = "none",
                 genesets = d$genesets, agg_method = "mean")
  fit <- run_engineering(X_train = X_na, params = params)

  expect_true(is.na(fit$X_transformed[1, "setA"]))
  expect_false(anyNA(fit$X_transformed[, "setB"]))
})

test_that("agg_method %in% c('pc1', 'ssgsea', 'gsva') error informatively on NA", {
  d <- .make_ssgsea_data()
  X_na <- d$X
  X_na[1, 1] <- NA

  for (m in c("pc1", "ssgsea", "gsva")) {
    params <- list(method = "engineer", col_transform = "none",
                   genesets = d$genesets, agg_method = m)
    expect_error(run_engineering(X_train = X_na, params = params), "NA")
  }
})
