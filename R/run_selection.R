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
#   - "variance"      : rank by feature variance (unsupervised)
#   - "pearson"       : rank by absolute univariate Pearson correlation with Y
#   - "spearman"      : rank by absolute univariate Spearman correlation with Y
#   - "relative_gain" : rank by CV prediction gain over a baseline model
#   - "rise"          : rank by adjusted p-value from RISE (SurrogateRank package)
# =============================================================================


# -----------------------------------------------------------------------------
#' Filter-based feature selection on a training matrix
#'
#' @description
#' Scores and ranks features (columns) of \code{X_train} according to a
#' filter criterion computed on the training data only, then retains either
#' the top \code{top_n} features or all features whose score meets
#' \code{threshold}. If both are supplied, \code{top_n} takes precedence.
#'
#' @details
#' Five filter methods are supported:
#' \itemize{
#'   \item \code{"variance"}: features are ranked by their variance across
#'     training samples. Does not use \code{Y_train}.
#'   \item \code{"pearson"}: features are ranked by the absolute value of
#'     their Pearson correlation with \code{Y_train}.
#'   \item \code{"spearman"}: features are ranked by the absolute value of
#'     their Spearman rank correlation with \code{Y_train}.
#'   \item \code{"relative_gain"}: features are ranked by their cross-validated
#'     prediction gain over a baseline model. For each feature j, an inner CV
#'     loop is run on \code{X_train}. In each inner fold, a baseline model
#'     (covariates only, or intercept-only if no covariates are provided) and
#'     a feature model (covariates + feature j) are fitted on the inner
#'     training partition and used to predict the inner test partition. The
#'     gain is the difference in predictive performance between the feature
#'     model and the baseline model, standardised so that positive gain always
#'     indicates that feature j improves prediction. Inner CV fold assignments
#'     are generated once and reused across all features, and the baseline
#'     predictions are computed once per inner fold. Supported metrics:
#'     \code{"rmse"}, \code{"srmse"}, \code{"r2"}, \code{"spearman"}.
#' }
#'
#' All scores are computed on \code{X_train} only. The selected feature names
#' are returned and used in \code{\link{predict_cv}} to subset both
#' \code{X_train} and \code{X_test} via
#' \code{X[, selected_features, drop = FALSE]}. No companion
#' \code{predict_selection} function is required.
#'
#' Scores for all features are returned in \code{selection_scores} for
#' diagnostic purposes.
#'
#' @param X_train Numeric matrix of dimensions n (samples) x p (features).
#'   Training predictor matrix. Must have column names.
#' @param Y_train Numeric vector of length n. Training response variable.
#'   Required for \code{"pearson"}, \code{"spearman"}, and
#'   \code{"relative_gain"}; ignored for \code{"variance"}.
#' @param covariates A numeric matrix of dimensions n x q to include in the
#'   baseline model for \code{"relative_gain"}. Ignored for all other methods.
#'   Pass \code{NULL} (default) for no covariates (intercept-only baseline).
#' @param treatment A binary numeric vector of length n with values 0 and 1,
#'   encoding treatment group membership. Required for \code{"rise"}; ignored
#'   for all other methods. Pass \code{NULL} (default) if not applicable.
#' @param params A named list of selection parameters with the following
#'   elements:
#'   \describe{
#'     \item{\code{method}}{Character string. One of \code{"variance"},
#'       \code{"pearson"}, \code{"spearman"}, or \code{"relative_gain"}.
#'       Required.}
#'     \item{\code{top_n}}{Positive integer. Number of top-ranked features to
#'       retain. Takes precedence over \code{threshold} if both are supplied.
#'       Either \code{top_n} or \code{threshold} must be specified.}
#'     \item{\code{threshold}}{Numeric. Minimum score a feature must achieve
#'       to be retained. Used only when \code{top_n} is \code{NULL}. For
#'       \code{"variance"}, a minimum variance; for \code{"pearson"} and
#'       \code{"spearman"}, a minimum absolute correlation; for
#'       \code{"relative_gain"}, a minimum gain; for \code{"rise"}, a
#'       maximum adjusted p-value (e.g. \code{0.05}).}
#'     \item{\code{rise_alpha}}{Numeric. Significance level passed to
#'       \code{rise.screen()} as \code{alpha}. Defaults to \code{0.05}.}
#'     \item{\code{rise_power_want_s}}{Numeric in (0,1). Desired power for
#'       surrogate test, passed as \code{power.want.s}. Either this or
#'       \code{rise_epsilon} must be specified.}
#'     \item{\code{rise_epsilon}}{Numeric in (0,1). Non-inferiority margin,
#'       passed as \code{epsilon}. Either this or \code{rise_power_want_s}
#'       must be specified.}
#'     \item{\code{rise_u_y_hyp}}{Numeric. Hypothesised treatment effect on
#'       the primary response on the probability scale, passed as
#'       \code{u.y.hyp}. Defaults to \code{NULL}.}
#'     \item{\code{rise_p_correction}}{Character. P-value adjustment method
#'       passed to \code{p.adjust()}, passed as \code{p.correction}.
#'       Defaults to \code{"BH"}.}
#'     \item{\code{rise_n_cores}}{Integer. Number of cores for parallel
#'       computation, passed as \code{n.cores}. Defaults to \code{1}.}
#'     \item{\code{rise_alternative}}{Character. Alternative hypothesis type,
#'       passed as \code{alternative}. One of \code{"less"} or
#'       \code{"two.sided"}. Defaults to \code{"two.sided"}.}
#'     \item{\code{rise_paired}}{Logical. Whether data are paired, passed as
#'       \code{paired}. Defaults to \code{FALSE}.}
#'     \item{\code{metric}}{Character string. Metric used to evaluate
#'       prediction quality in \code{"relative_gain"}. One of \code{"rmse"},
#'       \code{"srmse"}, \code{"r2"}, or \code{"spearman"}. Defaults to
#'       \code{"rmse"}. Ignored for other methods.}
#'     \item{\code{inner_folds}}{Positive integer. Number of inner CV folds
#'       for \code{"relative_gain"}. Defaults to \code{5}. Ignored for other
#'       methods.}
#'     \item{\code{seed}}{Integer. Random seed for inner fold assignment in
#'       \code{"relative_gain"}. Defaults to \code{12345}. Ignored for other
#'       methods.}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{selected_features}}{Character vector of selected column
#'       names, in decreasing order of score.}
#'     \item{\code{selection_method}}{Character string. The method used.}
#'     \item{\code{selection_scores}}{Named numeric vector of scores for
#'       \emph{all} features, in decreasing order. For
#'       \code{"relative_gain"}, scores are the gain values.}
#'     \item{\code{top_n}}{Integer or \code{NULL}. The \code{top_n} value
#'       used after resolving precedence with \code{threshold}.}
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
#'
#' # Relative gain: retain features with positive gain over null model
#' result <- run_selection(X, Y,
#'   params = list(method = "relative_gain", threshold = 0,
#'                 metric = "rmse", inner_folds = 5))
#'
#' # Relative gain with covariates
#' covariates <- matrix(rnorm(40 * 2), nrow = 40,
#'                      dimnames = list(NULL, c("age", "sex")))
#' result <- run_selection(X, Y, covariates = covariates,
#'   params = list(method = "relative_gain", threshold = 0, metric = "r2"))
#' }
#'
#' @export
# -----------------------------------------------------------------------------
run_selection <- function(X_train, Y_train = NULL, covariates = NULL,
                          treatment = NULL, params) {

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
  if (method %in% c("pearson", "spearman", "relative_gain", "rise")) {
    if (is.null(Y_train))
      stop("[predictomics] Y_train must be provided for method = '",
           method, "'.", call. = FALSE)
    .validate_Y(Y_train)
    .validate_Y_X_compat(Y_train, X_train)
  }

  if (method == "rise") {
    if (is.null(treatment))
      stop("[predictomics] treatment must be provided for method = 'rise'.",
           call. = FALSE)
    if (!is.numeric(treatment) || !all(treatment %in% c(0, 1)))
      stop("[predictomics] treatment must be a binary numeric vector (0/1) ",
           "for method = 'rise'.", call. = FALSE)
    if (!requireNamespace("SurrogateRank", quietly = TRUE))
      stop("[predictomics] The SurrogateRank package is required for ",
           "method = 'rise'. Install it with: install.packages('SurrogateRank')",
           call. = FALSE)
  }

  # top_n takes precedence - inform the user
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
                   },

                   relative_gain = {
                     .compute_relative_gain(
                       X_train      = X_train,
                       Y_train      = Y_train,
                       covariates   = covariates,
                       metric       = params$metric       %||% "rmse",
                       inner_folds  = params$inner_folds  %||% 5L,
                       seed         = params$seed         %||% 12345L
                     )
                   },

                   rise = {
                     .compute_rise_scores(
                       X_train          = X_train,
                       Y_train          = Y_train,
                       treatment        = treatment,
                       top_n            = top_n,
                       alpha            = params$rise_alpha         %||% 0.05,
                       power_want_s     = params$rise_power_want_s,
                       epsilon          = params$rise_epsilon,
                       u_y_hyp          = params$rise_u_y_hyp       %||% NULL,
                       p_correction     = params$rise_p_correction  %||% "BH",
                       n_cores          = params$rise_n_cores        %||% 1L,
                       alternative      = params$rise_alternative   %||% "two.sided",
                       paired           = params$rise_paired         %||% FALSE
                     )
                   }
  )

  # Sort descending for reporting and selection
  scores <- sort(scores, decreasing = TRUE)

  # ---------------------------------------------------------------------------
  # 3. Select features
  # ---------------------------------------------------------------------------
  selected_features <- if (!is.null(top_n)) {

    if (method == "relative_gain" & any(scores[seq_len(top_n)] < 0)){

      n_below_floor = sum(scores[seq_len(top_n)] < 0)

      message(
        "[predictomics] Relative gain selection: ", n_below_floor,
        " feature(s) in the top ", top_n,
        " have relative gain < 0.",
        " By definition, these features negatively affect the predictive performance",
        " and have therefore been removed from selection."
      )

    names(scores)[seq_len(top_n)][which(scores[seq_len(top_n)] > 0)]

    } else {

    names(scores)[seq_len(top_n)]

    }

  } else {

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
#' Compute univariate CV relative gain scores for all features
#'
#' @description
#' For each feature j in \code{X_train}, runs an inner CV loop to compute the
#' predictive gain of adding feature j to a baseline model. The baseline model
#' contains \code{covariates} only (or an intercept if \code{covariates} is
#' \code{NULL}). Inner fold assignments are generated once and shared across
#' all features. Baseline predictions are computed once per inner fold.
#'
#' Gain is defined so that positive values always indicate that feature j
#' improves prediction over baseline:
#' \itemize{
#'   \item \code{"rmse"}, \code{"srmse"}: gain = baseline metric - feature
#'     metric (lower is better, so positive gain = improvement).
#'   \item \code{"r2"}, \code{"spearman"}: gain = feature metric - baseline
#'     metric (higher is better, so positive gain = improvement).
#' }
#'
#' @param X_train Numeric matrix. Training features.
#' @param Y_train Numeric vector. Training response.
#' @param covariates Numeric matrix or \code{NULL}. Baseline covariates.
#' @param metric Character string. One of \code{"rmse"}, \code{"srmse"},
#'   \code{"r2"}, \code{"spearman"}.
#' @param inner_folds Integer. Number of inner CV folds.
#' @param seed Integer. Seed for inner fold assignment.
#'
#' @return Named numeric vector of gain scores, one per feature.
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_relative_gain <- function(X_train, Y_train, covariates, metric,
                                   inner_folds, seed) {

  n         <- nrow(X_train)
  p         <- ncol(X_train)
  feat_names <- colnames(X_train)

  # ---------------------------------------------------------------------------
  # Generate inner fold assignments once, shared across all features
  # ---------------------------------------------------------------------------
  inner_fold_ids <- make_folds(n = n, cv_type = "kfold",
                               k = inner_folds, seed = seed)

  # ---------------------------------------------------------------------------
  # Pre-build baseline design matrix (covariates or intercept only)
  # Used identically across all feature models
  # ---------------------------------------------------------------------------
  has_covariates <- !is.null(covariates) && ncol(covariates) > 0L

  # ---------------------------------------------------------------------------
  # Compute baseline CV predictions once (shared across all features)
  # ---------------------------------------------------------------------------
  baseline_pred <- numeric(n)

  for (f in seq_len(inner_folds)) {

    tr  <- which(inner_fold_ids != f)
    tst <- which(inner_fold_ids == f)

    Y_tr  <- Y_train[tr]
    Y_tst <- Y_train[tst]

    if (has_covariates) {
      cov_tr  <- as.data.frame(covariates[tr,  , drop = FALSE])
      cov_tst <- as.data.frame(covariates[tst, , drop = FALSE])
      df_tr   <- cbind(data.frame(.Y = Y_tr), cov_tr)
      df_tst  <- cov_tst
      fit     <- lm(.Y ~ ., data = df_tr)
    } else {
      df_tr  <- data.frame(.Y = Y_tr)
      df_tst <- data.frame(.intercept = rep(1, length(tst)))
      fit    <- lm(.Y ~ 1, data = df_tr)
    }

    baseline_pred[tst] <- predict(fit, newdata = df_tst)
  }

  baseline_score <- .compute_metric(Y_train, baseline_pred, metric)

  # ---------------------------------------------------------------------------
  # Compute feature model CV predictions and gain for each feature
  # ---------------------------------------------------------------------------
  gains <- numeric(p)
  names(gains) <- feat_names

  for (j in seq_len(p)) {

    feat_pred <- numeric(n)

    for (f in seq_len(inner_folds)) {

      tr  <- which(inner_fold_ids != f)
      tst <- which(inner_fold_ids == f)

      Y_tr   <- Y_train[tr]
      feat_j <- X_train[, j]

      if (has_covariates) {
        cov_tr  <- as.data.frame(covariates[tr,  , drop = FALSE])
        cov_tst <- as.data.frame(covariates[tst, , drop = FALSE])
        df_tr   <- cbind(data.frame(.Y = Y_tr, .feat = feat_j[tr]),  cov_tr)
        df_tst  <- cbind(data.frame(.feat = feat_j[tst]),             cov_tst)
      } else {
        df_tr  <- data.frame(.Y = Y_tr,       .feat = feat_j[tr])
        df_tst <- data.frame(.feat = feat_j[tst])
      }

      fit <- lm(.Y ~ ., data = df_tr)
      feat_pred[tst] <- predict(fit, newdata = df_tst)
    }

    feature_score <- .compute_metric(Y_train, feat_pred, metric)
    gains[j]      <- .compute_gain(baseline_score, feature_score, metric)
  }

  gains
}


# -----------------------------------------------------------------------------
#' Compute a scalar prediction metric from observed and predicted vectors
#'
#' @param obs Numeric vector of observed values.
#' @param pred Numeric vector of predicted values.
#' @param metric Character string. One of \code{"rmse"}, \code{"srmse"},
#'   \code{"r2"}, \code{"spearman"}.
#' @return A single numeric value.
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_metric <- function(obs, pred, metric) {

  switch(metric,
         rmse     = sqrt(mean((obs - pred)^2)),
         srmse    = sqrt(mean((obs - pred)^2)) / sd(obs),
         r2       = cor(obs, pred, method = "pearson")^2,
         spearman = cor(obs, pred, method = "spearman")
  )
}


# -----------------------------------------------------------------------------
#' Compute directional gain between baseline and feature model scores
#'
#' @description
#' Returns gain such that positive values always indicate improvement of the
#' feature model over baseline, regardless of metric direction.
#'
#' @param baseline_score Numeric. Metric value for the baseline model.
#' @param feature_score Numeric. Metric value for the feature model.
#' @param metric Character string. The metric used.
#' @return A single numeric gain value.
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_gain <- function(baseline_score, feature_score, metric) {

  # For lower-is-better metrics: gain = baseline - feature (positive = better)
  # For higher-is-better metrics: gain = feature - baseline (positive = better)
  if (metric %in% c("rmse", "srmse")) {
    baseline_score - feature_score
  } else {
    feature_score - baseline_score
  }
}


# -----------------------------------------------------------------------------
#' Compute RISE adjusted p-value scores for feature selection
#'
#' @description
#' Calls \code{SurrogateRank::rise.screen()} on the training data, reshaping
#' inputs from the predictomics format (unified \code{X_train}, \code{Y_train},
#' binary \code{treatment}) into the RISE format (\code{yone}, \code{yzero},
#' \code{sone}, \code{szero}). Returns a named numeric vector of negated
#' adjusted p-values (so that sorting descending yields features with the
#' smallest p-values first, consistent with the rest of \code{run_selection}).
#'
#' When \code{top_n} is non-NULL and more than \code{top_n} features have
#' adjusted p-values of exactly 1 (a common occurrence due to multiplicity
#' penalties in high-dimensional settings), a note is printed and unadjusted
#' p-values are used as tiebreakers among features with adjusted p-value = 1.
#'
#' @param X_train Numeric matrix. Training predictor matrix.
#' @param Y_train Numeric vector. Training response.
#' @param treatment Binary numeric vector (0/1). Treatment assignment.
#' @param top_n Integer or NULL. Used only to determine whether the tiebreak
#'   note should be printed.
#' @param alpha,power_want_s,epsilon,u_y_hyp,p_correction,n_cores,alternative,paired
#'   Arguments passed directly to \code{SurrogateRank::rise.screen()}.
#'
#' @return Named numeric vector of negated adjusted p-values, one per feature.
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_rise_scores <- function(X_train, Y_train, treatment, top_n,
                                 alpha, power_want_s, epsilon, u_y_hyp,
                                 p_correction, n_cores, alternative, paired) {

  # Emit pairing order note when paired = TRUE
  if (isTRUE(paired)) {
    message(
      "[predictomics] RISE paired mode: assuming samples are in matched ",
      "pre/post order - row i of the pre-treatment group (treatment == 0) ",
      "corresponds to row i of the post-treatment group (treatment == 1). ",
      "Ensure this ordering is correct before proceeding."
    )
  }

  # ---------------------------------------------------------------------------
  # 1. Reshape inputs into RISE format
  # ---------------------------------------------------------------------------
  idx1  <- which(treatment == 1)
  idx0  <- which(treatment == 0)

  yone  <- Y_train[idx1]
  yzero <- Y_train[idx0]
  sone  <- X_train[idx1, , drop = FALSE]
  szero <- X_train[idx0, , drop = FALSE]

  # ---------------------------------------------------------------------------
  # 2. Call rise.screen(), suppressing its internal plot and verbose output
  # ---------------------------------------------------------------------------
  res <- SurrogateRank::rise.screen(
    yone               = yone,
    yzero              = yzero,
    sone               = sone,
    szero              = szero,
    alpha              = alpha,
    power.want.s       = power_want_s,
    epsilon            = epsilon,
    u.y.hyp            = u_y_hyp,
    p.correction       = p_correction,
    n.cores            = n_cores,
    alternative        = alternative,
    paired             = paired,
    return.all.screen  = TRUE,
    return.screen.plot = FALSE,
    return.all.weights = FALSE,
    verbose            = FALSE
  )

  metrics <- res[["screening.metrics"]]

  # Align to X_train column order
  p_adj   <- setNames(metrics$p_adjusted,   metrics$marker)
  p_unadj <- setNames(metrics$p_unadjusted, metrics$marker)
  p_adj   <- p_adj[colnames(X_train)]
  p_unadj <- p_unadj[colnames(X_train)]

  # ---------------------------------------------------------------------------
  # 3. Handle ceiling at 1: use unadjusted p-values as tiebreaker
  # ---------------------------------------------------------------------------
  n_below_ceiling <- sum(p_adj < 1)

  if (!is.null(top_n) && n_below_ceiling < top_n) {
    message(
      "[predictomics] RISE: ", n_below_ceiling, " feature(s) have adjusted ",
      "p-values < 1 (out of ", length(p_adj), " total). Since top_n = ", top_n,
      " exceeds this, unadjusted p-values will be used as tiebreakers for ",
      "features with adjusted p-value = 1. Rankings among these tiebroken ",
      "features are not meaningful in the adjusted sense."
    )
  }

  # ---------------------------------------------------------------------------
  # 4. Construct composite score: negate so that sort(decreasing=TRUE) gives
  #    lowest-p features first. Features with adj p < 1 are ranked by -p_adj;
  #    features at the ceiling are ranked by -p_unadj, offset to be strictly
  #    below all non-ceiling scores.
  # ---------------------------------------------------------------------------
  ceiling_offset  <- if (any(p_adj < 1)) min(-p_adj[p_adj < 1]) - 1 else -1

  scores <- ifelse(
    p_adj < 1,
    -p_adj,
    ceiling_offset - p_unadj   # pushes ceiling features below all non-ceiling
  )

  names(scores) <- colnames(X_train)
  scores
}
