# -----------------------------------------------------------------------------
# plot.predictomics() highlighting
# -----------------------------------------------------------------------------

.make_plot_data <- function(n = 30, p = 10, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("gene", seq_len(p))
  Y <- X[, 1] * 2 + rnorm(n)
  list(X = X, Y = Y)
}

test_that("plot.predictomics works without highlight or treatment", {
  d <- .make_plot_data()
  result <- predict_cv(Y = d$Y, X = d$X, folds = 3, verbose = FALSE)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})

test_that("highlight requires highlight_type", {
  d <- .make_plot_data()
  result <- predict_cv(Y = d$Y, X = d$X, folds = 3, verbose = FALSE)

  expect_error(
    plot(result, highlight = rep(c("A", "B"), length.out = length(result$observed))),
    "highlight_type"
  )

  expect_error(
    plot(result, highlight = rep(c("A", "B"), length.out = length(result$observed)),
        highlight_type = "not_a_type"),
    "highlight_type"
  )
})

test_that("highlight must match the length of x$observed", {
  d <- .make_plot_data()
  result <- predict_cv(Y = d$Y, X = d$X, folds = 3, verbose = FALSE)

  expect_error(
    plot(result, highlight = c("A", "B"), highlight_type = "categorical"),
    "length"
  )
})

test_that("categorical highlight with <= 8 levels uses a Set2 legend", {
  d <- .make_plot_data()
  result <- predict_cv(Y = d$Y, X = d$X, folds = 3, verbose = FALSE)
  n <- length(result$observed)
  batch <- factor(rep(c("A", "B", "C"), length.out = n))

  p <- plot(result, highlight = batch, highlight_type = "categorical",
           highlight_label = "Batch")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$colour, "Batch")
  expect_true(any(vapply(p$scales$scales, function(s) "colour" %in% s$aesthetics,
                        logical(1))))
})

test_that("categorical highlight with > 8 levels uses a discrete viridis scale", {
  d <- .make_plot_data(n = 40)
  result <- predict_cv(Y = d$Y, X = d$X, folds = 3, verbose = FALSE)
  n <- length(result$observed)
  many_levels <- factor(paste0("g", seq_len(n) %% 10))

  p <- plot(result, highlight = many_levels, highlight_type = "categorical")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$colour, "Highlight")
})

test_that("continuous highlight uses a viridis gradient with a colourbar legend", {
  d <- .make_plot_data()
  result <- predict_cv(Y = d$Y, X = d$X, folds = 3, verbose = FALSE)
  n <- length(result$observed)
  age <- rnorm(n)

  p <- plot(result, highlight = age, highlight_type = "continuous",
           highlight_label = "Age")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$colour, "Age")
  expect_true(is.numeric(p$data$highlight_value))
})

test_that("highlight overrides treatment colouring entirely", {
  d <- .make_plot_data()
  treatment <- rep(c(0, 1), length.out = length(d$Y))
  result <- predict_cv(Y = d$Y, X = d$X, treatment = treatment, folds = 3,
                       verbose = FALSE)
  n <- length(result$observed)
  batch <- factor(rep(c("A", "B"), length.out = n))

  p <- plot(result, highlight = batch, highlight_type = "categorical")

  expect_equal(p$labels$colour, "Highlight")
  expect_false("treatment_group" %in% names(p$data))
  expect_true("highlight_value" %in% names(p$data))
})
