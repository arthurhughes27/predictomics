# =============================================================================
# metrics.R
# Performance metrics for predictomics result objects.
# =============================================================================


# -----------------------------------------------------------------------------
#' Compute prediction performance metrics
#'
#' @description
#' Computes prediction performance metrics from a \code{predictomics} result
#' object returned by \code{\link{predict_cv}}. Metrics are computed by
#' comparing cross-validated predicted values to observed values.
#'
#' @details
#' The following metrics are computed:
#' \itemize{
#'   \item \strong{RMSE}: Root mean squared error,
#'     \eqn{\sqrt{\frac{1}{n}\sum_{i=1}^{n}(y_i - \hat{y}_i)^2}}.
#'   \item \strong{sRMSE}: Standardised RMSE, defined as RMSE divided by the
#'     standard deviation of the observed values. Values below 1 indicate that
#'     the model outperforms a naive mean predictor.
#'   \item \strong{R2}: Squared Pearson correlation between observed and
#'     predicted values, \eqn{r^2_{pearson}}. Bounded between 0 and 1.
#'     Note that this differs from the coefficient of determination
#'     (1 - SS_res / SS_tot), which can be negative.
#'   \item \strong{SpearmanR}: Spearman rank correlation between observed and
#'     predicted values. Robust to outliers and monotone nonlinear
#'     relationships.
#' }
#'
#' @param x A \code{predictomics} object returned by \code{\link{predict_cv}}.
#' @param digits Integer. Number of decimal places to round metrics to.
#'   Defaults to \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector containing \code{RMSE}, \code{sRMSE},
#'   \code{R2}, and \code{SpearmanR}.
#'
#' @seealso \code{\link{predict_cv}}, \code{\link{plot.predictomics}},
#'   \code{\link{print.predictomics}}
#'
#' @examples
#' \dontrun{
#' result <- predict_cv(Y = Y, X = X, model_params = list(method = "glmnet"))
#' metrics(result)
#' }
#'
#' @export
# -----------------------------------------------------------------------------
metrics <- function(x, ...) UseMethod("metrics")


# -----------------------------------------------------------------------------
#' @rdname metrics
#' @export
# -----------------------------------------------------------------------------
metrics.predictomics <- function(x, digits = 4, ...) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!inherits(x, "predictomics"))
    stop("[predictomics] x must be a predictomics object returned by ",
         "predict_cv().", call. = FALSE)

  if (!is.numeric(digits) || length(digits) != 1L || digits < 0L)
    stop("[predictomics] digits must be a non-negative integer.", call. = FALSE)

  obs  <- x$observed
  pred <- x$predicted
  n    <- length(obs)

  # ---------------------------------------------------------------------------
  # 2. Compute metrics
  # ---------------------------------------------------------------------------
  rmse     <- sqrt(mean((obs - pred)^2))
  srmse    <- rmse / sd(obs)
  r2       <- cor(obs, pred, method = "pearson")^2
  spearman <- cor(obs, pred, method = "spearman")

  # ---------------------------------------------------------------------------
  # 3. Return
  # ---------------------------------------------------------------------------
  round(
    c(RMSE = rmse, sRMSE = srmse, R2 = r2, SpearmanR = spearman),
    digits = digits
  )
}
