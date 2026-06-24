# =============================================================================
# run_selection.R
# Feature selection for the predictomics package.
#
# Implements filter-based variable selection methods that score and rank
# features based on statistics computed on the training set only. Selected
# feature names are returned for direct application to the test set via
# column subsetting in predict_cv().
#
# Supported methods:
#   - "variance"    : rank by feature variance (unsupervised)
#   - "pearson"     : rank by absolute univariate Pearson correlation with Y
#   - "spearman"    : rank by absolute univariate Spearman correlation with Y
# =============================================================================


# -----------------------------------------------------------------------------
#' Filter-based feature selection on a training matrix
#'
#' @description
#' Scores and ranks features (columns) of \code{X_train} according to a
#' filter criterion computed on the training data only, then retains either
#' the top \code{top_n} features or all features whose score exceeds
#' \code{threshold}. If both are supplied, \code{top_n} takes precedence.
#'
#' @details
#' Three filter methods are supported:
#' \itemize{
#'   \item \code{"variance"}: features are ranked by their variance across
#'     training samples, computed via \code{var()}. Does not use \code{Y_train}.
#'   \item \code{"pearson"}: features are ranked by the absolute value of their
#'     Pearson correlation with \code{Y_train}.
#'   \item \code{"spearman"}: features are ranked by the absolute value of their
#'     Spearman rank correlation with \code{Y_train}.
#' }
#'
#' All scores are computed on \code{X_train} only. The selected feature names
#' are returned and used in \code{\link{predict_cv}} to subset both
#' \code{X_train} and \code{X_test} via
#' \code{X[, selected_features, drop = FALSE]}. No companion
#' \code{predict_selection} function is required.
#'
#' Scores for all features (not only selected ones) are returned in
#' \code{selection_scores} for diagnostic purposes.
#'
#' @param X_train Numeric matrix of dimensions n (samples) x p (features).
#'   Training predictor matrix. Must have column names.
#' @param Y_train Numeric vector of length n. Training response variable.
#'   Required for \code{"pearson"} and \code{"spearman"}; ignored for
#'   \code{"variance"}.
#' @param covariates A numeric matrix of dimensions n (samples) x q
#'   (covariates) to include in the baseline model for supervised filter
#'   methods that require a baseline (e.g. \code{"relative_gain"}). Ignored
#'   for \code{"variance"}, \code{"pearson"}, and \code{"spearman"}.
#'   Pass \code{NULL} (default) for no covariates.
#' @param params A named list of selection parameters with the following
#'   elements:
#'   \describe{
#'     \item{\code{method}}{Character string. One of \code{"variance"},
#'       \code{"pearson"}, or \code{"spearman"}. Required.}
#'     \item{\code{top_n}}{Positive integer. Number of top-ranked features to
#'       retain. Takes precedence over \code{threshold} if both are supplied.
#'       Either \code{top_n} or \code{threshold} must be specified.}
#'     \item{\code{threshold}}{Numeric. Minimum score a feature must achieve
#'       to be retained. Used only when \code{top_n} is \code{NULL}. For
#'       \code{"variance"}, this is a minimum variance; for \code{"pearson"}
#'       and \code{"spearman"}, a minimum absolute correlation.}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{selected_features}}{Character vector of selected column
#'       names, in decreasing order of score.}
#'     \item{\code{selection_method}}{Character string. The method used.}
#'     \item{\code{selection_scores}}{Named numeric vector of scores for
#'       \emph{all} features, in decreasing order. Useful for diagnostics
#'       and inspection.}
#'     \item{\code{top_n}}{Integer or \code{NULL}. The \code{top_n} value
#'       used, after resolving precedence with \code{threshold}.}
#'     \item{\code{threshold}}{Numeric or \code{NULL}. The \code{threshold}
#'       value used, or \code{NULL} if \code{top_n} was applied.}
#'   }
#'
#' @seealso \code{\link{predict_cv}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(40 * 200), nrow = 40)
#' colnames(X) <- paste0("gene", 1:200)
#' Y <- X[, 1] * 2 + rnorm(40)
#'
#' # Variance filter: retain top 50 features
#' result <- run_selection(X, Y,
#'   params = list(method = "variance", top_n = 50))
#' length(result$selected_features)  # 50
#'
#' # Pearson correlation filter: retain features with |r| > 0.2
#' result <- run_selection(X, Y,
#'   params = list(method = "pearson", threshold = 0.2))
#'
#' # Spearman correlation filter: top 100 features
#' # (top_n takes precedence over threshold when both supplied)
#' result <- run_selection(X, Y,
#'   params = list(method = "spearman", top_n = 100, threshold = 0.1))
#' }
#'
#' @export
# -----------------------------------------------------------------------------
run_selection <- function(X_train, Y_train = NULL, covariates = NULL, params) {

  # ---------------------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------------------
  .validate_X(X_train)
  .validate_selection_params(params, p = ncol(X_train))

  if (is.null(colnames(X_train)))
    stop("[predictomics] X_train must have column names for feature selection.",
         call. = FALSE)

  method    <- params$method
  top_n     <- params$top_n
  threshold <- params$threshold

  # Supervised methods require Y_train
  if (method %in% c("pearson", "spearman")) {
    if (is.null(Y_train))
      stop("[predictomics] Y_train must be provided for method = '", method,
           "'.", call. = FALSE)
    .validate_Y(Y_train)
    .validate_Y_X_compat(Y_train, X_train)
  }

  # top_n takes precedence — log a message so the user is aware
  if (!is.null(top_n) && !is.null(threshold)) {
    message("[predictomics] Both top_n and threshold supplied to run_selection; ",
            "top_n takes precedence.")
    threshold <- NULL
  }

  # ---------------------------------------------------------------------------
  # 2. Compute scores
  # ---------------------------------------------------------------------------
  scores <- switch(method,

                   variance = {
                     apply(X_train, 2, var)
                   },

                   pearson = {
                     apply(X_train, 2, function(x) abs(cor(x, Y_train, method = "pearson")))
                   },

                   spearman = {
                     apply(X_train, 2, function(x) abs(cor(x, Y_train, method = "spearman")))
                   }
  )

  # Sort descending for reporting and selection
  scores <- sort(scores, decreasing = TRUE)

  # ---------------------------------------------------------------------------
  # 3. Select features
  # ---------------------------------------------------------------------------
  selected_features <- if (!is.null(top_n)) {

    names(scores)[seq_len(top_n)]

  } else {

    # threshold mode
    sel <- names(scores)[scores >= threshold]

    if (length(sel) == 0L)
      stop("[predictomics] No features pass the threshold (", threshold,
           ") for method = '", method, "'. ",
           "Consider lowering threshold or using top_n instead.",
           call. = FALSE)

    sel
  }

  # ---------------------------------------------------------------------------
  # 4. Return
  # ---------------------------------------------------------------------------
  list(
    selected_features = selected_features,
    selection_method  = method,
    selection_scores  = scores,
    top_n             = top_n,
    threshold         = threshold
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# -----------------------------------------------------------------------------
#' Validate selection_params
#'
#' @description
#' Checks that \code{selection_params} specifies a supported method, that
#' exactly one of \code{top_n} or \code{threshold} is specified (or both,
#' since \code{top_n} takes precedence), and that their values are valid
#' relative to the number of available features \code{p}.
#'
#' @param params The \code{selection_params} list.
#' @param p Integer. Number of features in \code{X_train}, used to bound
#'   \code{top_n}.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_selection_params <- function(params, p) {

  # method
  supported <- c("variance", "pearson", "spearman")
  if (!params$method %in% supported)
    stop("[predictomics] selection_params$method must be one of: ",
         paste(supported, collapse = ", "), ".", call. = FALSE)

  top_n     <- params$top_n
  threshold <- params$threshold

  # At least one of top_n or threshold must be provided
  if (is.null(top_n) && is.null(threshold))
    stop("[predictomics] selection_params must specify at least one of ",
         "'top_n' or 'threshold'.", call. = FALSE)

  # Validate top_n if provided
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L ||
        top_n != as.integer(top_n) || top_n < 1L)
      stop("[predictomics] selection_params$top_n must be a positive integer.",
           call. = FALSE)
    if (top_n > p)
      stop("[predictomics] selection_params$top_n (", top_n, ") exceeds the ",
           "number of available features (", p, ").", call. = FALSE)
  }

  # Validate threshold if provided
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0)
      stop("[predictomics] selection_params$threshold must be a single ",
           "non-negative numeric value.", call. = FALSE)
  }

  invisible(NULL)
}
