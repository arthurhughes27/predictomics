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
#   - "dearseq"       : rank by adjusted p-value from a dearseq differential
#                       expression analysis (dearseq package)
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
#' If \code{top_n} is greater than or equal to the number of available
#' features (\code{ncol(X_train)}, or the number of \code{genesets} for
#' \code{method = "dearseq"} with \code{dearseq_level = "geneset"}), every
#' feature (or geneset) is selected and a message is printed - selection
#' scores are \strong{not} computed in this case, since every candidate would
#' be selected regardless of its score. This is checked before any scoring
#' computation begins, so no method-specific computation (e.g. inner-CV
#' folds for \code{"relative_gain"}, or an external call for \code{"rise"}/
#' \code{"dearseq"}) is performed when it would have no effect on the
#' selection outcome. \code{selection_scores} is a vector of \code{NA} in
#' this case.
#'
#' Five filter methods are supported:
#' \itemize{
#'   \item \code{"variance"}: features are ranked by their variance across
#'     training samples. Does not use \code{Y_train}. When called from
#'     \code{\link{predict_cv}} with \code{engineering_params$col_transform =
#'     "z"}, \code{X_train} here is automatically a separately re-engineered,
#'     pre-z-score matrix rather than the fully engineered one - see
#'     \code{\link{predict_cv}}'s Details for why.
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
#'   \item \code{"dearseq"}: features are ranked by adjusted p-value from a
#'     dearseq differential expression analysis (requires the Suggested
#'     \pkg{dearseq} package). Two independent choices control the analysis:
#'     \itemize{
#'       \item \code{dearseq_mode} chooses \emph{which comparison} is run:
#'         \code{"classic"} (default) compares \code{treatment == 1} vs
#'         \code{treatment == 0} (\code{treatment} required); if
#'         \code{individual_id} and \code{timepoint} are also supplied,
#'         \code{X_train} is first restricted to \code{timepoint == 1} rows
#'         before the comparison is run. \code{"paired"} compares
#'         \code{timepoint == 1} vs \code{timepoint == 0}
#'         (\code{individual_id} and \code{timepoint} required, with every
#'         individual having exactly one observation at each timepoint),
#'         with \code{sample_group = individual_id}; if \code{treatment} is
#'         also supplied, \code{X_train} is first restricted to
#'         \code{treatment == 1} rows before the comparison is run.
#'       \item \code{dearseq_level} chooses \emph{what is tested}:
#'         \code{"gene"} (default) calls \code{dearseq::dear_seq()} and
#'         scores every column of \code{X_train} individually. \code{"geneset"}
#'         calls \code{dearseq::dgsa_seq()} with \code{genesets} instead, and
#'         scores each geneset as a unit; the selected genesets (by
#'         \code{top_n}/\code{threshold} on their adjusted p-values) are then
#'         resolved to the union of their member genes (restricted to
#'         columns of \code{X_train}) for \code{selected_features}, while
#'         \code{selection_scores} remains indexed by geneset name.
#'     }
#'     In all cases, only \emph{rows} are ever restricted before the
#'     comparison; adjusted p-values are computed from the (possibly
#'     row-restricted) comparison but scoring always considers every gene
#'     (or geneset) supplied.
#' }
#'
#' **Missing values.** \code{"variance"}, \code{"pearson"}, and
#' \code{"spearman"} tolerate \code{NA} values in \code{X_train} (variance is
#' computed with \code{na.rm = TRUE}; correlations use
#' \code{use = "pairwise.complete.obs"}). \code{"relative_gain"},
#' \code{"rise"}, and \code{"dearseq"} do not support \code{NA} and will
#' error if any is present in \code{X_train}.
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
#'   baseline model for \code{"relative_gain"}. Also used as the default
#'   covariates for \code{"dearseq"}, unless \code{params$dearseq_covariates}
#'   is supplied. Ignored for all other methods. Pass \code{NULL} (default)
#'   for no covariates (intercept-only baseline).
#' @param treatment A binary numeric vector of length n with values 0 and 1,
#'   encoding treatment group membership. Required for \code{"rise"} (unless
#'   \code{rise_paired = TRUE}, in which case \code{individual_id}/
#'   \code{timepoint} are required instead) and for \code{"dearseq"} with
#'   \code{dearseq_mode = "classic"}; optionally used to restrict rows for
#'   \code{"dearseq"} with \code{dearseq_mode = "paired"}. Ignored for all
#'   other methods. Pass \code{NULL} (default) if not applicable.
#' @param individual_id Vector of length n identifying individuals. Required
#'   for \code{"rise"} with \code{rise_paired = TRUE} and for \code{"dearseq"}
#'   with \code{dearseq_mode = "paired"} (used as \code{sample_group});
#'   optionally used (together with \code{timepoint}) to restrict rows for
#'   \code{dearseq_mode = "classic"}. Ignored for all other methods. Pass
#'   \code{NULL} (default) if not applicable.
#' @param timepoint A binary numeric vector of length n (0/1), paired with
#'   \code{individual_id}. Required for \code{"rise"} with
#'   \code{rise_paired = TRUE} (used as the pre/post contrast in place of
#'   \code{treatment}) and for \code{"dearseq"} with
#'   \code{dearseq_mode = "paired"}; optionally used to restrict rows for
#'   \code{dearseq_mode = "classic"}. Ignored for all other methods. Pass
#'   \code{NULL} (default) if not applicable.
#' @param params A named list of selection parameters with the following
#'   elements:
#'   \describe{
#'     \item{\code{method}}{Character string. One of \code{"variance"},
#'       \code{"pearson"}, \code{"spearman"}, \code{"relative_gain"},
#'       \code{"rise"}, or \code{"dearseq"}. Required.}
#'     \item{\code{top_n}}{Positive integer. Number of top-ranked features to
#'       retain. Takes precedence over \code{threshold} if both are supplied.
#'       Either \code{top_n} or \code{threshold} must be specified. If
#'       \code{top_n} is greater than or equal to the number of available
#'       features (or genesets, for \code{dearseq_level = "geneset"}), all of
#'       them are selected and a message is printed; selection scores are not
#'       computed in this case (see Details).}
#'     \item{\code{threshold}}{Numeric. Minimum score a feature must achieve
#'       to be retained. Used only when \code{top_n} is \code{NULL}. For
#'       \code{"variance"}, a minimum variance; for \code{"pearson"} and
#'       \code{"spearman"}, a minimum absolute correlation; for
#'       \code{"relative_gain"}, a minimum gain; for \code{"rise"} and
#'       \code{"dearseq"}, a maximum adjusted p-value (e.g. \code{0.05}).}
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
#'       \code{paired}. When \code{TRUE}, \code{individual_id} and
#'       \code{timepoint} are required (in the same paired design as
#'       \code{"dearseq"} with \code{dearseq_mode = "paired"}): the
#'       pre/post contrast is \code{timepoint} rather than \code{treatment},
#'       and rows are sorted by \code{individual_id} (then \code{timepoint})
#'       to guarantee correct pairing regardless of input row order.
#'       Defaults to \code{FALSE}.}
#'     \item{\code{dearseq_mode}}{Character string. One of \code{"classic"}
#'       (default) or \code{"paired"}. See Details.}
#'     \item{\code{dearseq_level}}{Character string. One of \code{"gene"}
#'       (default) or \code{"geneset"}. See Details.}
#'     \item{\code{genesets}}{Named list of character vectors of feature
#'       names, in the same form as \code{engineering_params$genesets}.
#'       Required when \code{dearseq_level = "geneset"}; ignored otherwise.}
#'     \item{\code{dearseq_which_test}}{Character string. Passed to
#'       \code{dear_seq()} as \code{which_test}. One of \code{"asymptotic"}
#'       (default) or \code{"permutation"}.}
#'     \item{\code{dearseq_preprocessed}}{Logical. Passed to \code{dear_seq()}
#'       as \code{preprocessed}. Defaults to \code{TRUE}.}
#'     \item{\code{dearseq_padjust_methods}}{Character string. Passed to
#'       \code{dear_seq()} as \code{padjust_methods}. One of
#'       \code{stats::p.adjust.methods}. Defaults to \code{"BH"}.}
#'     \item{\code{dearseq_which_weights}}{Character string. Passed to
#'       \code{dear_seq()} as \code{which_weights}. One of \code{"loclin"}
#'       (default), \code{"voom"}, or \code{"none"}.}
#'     \item{\code{dearseq_n_perm}}{Positive integer. Passed to
#'       \code{dear_seq()} as \code{n_perm}. Defaults to \code{1000}.}
#'     \item{\code{dearseq_bw}}{Character string or positive numeric. Passed
#'       to \code{dear_seq()} as \code{bw}. Defaults to \code{"nrd"}.}
#'     \item{\code{dearseq_kernel}}{Character string. Passed to
#'       \code{dear_seq()} as \code{kernel}. Defaults to \code{"gaussian"}.}
#'     \item{\code{dearseq_covariates}}{A numeric matrix or data frame of
#'       covariates to use specifically for the dearseq comparison, in place
#'       of the \code{covariates} argument. Useful when the covariates
#'       desired for differential expression testing differ from those used
#'       in the main prediction pipeline. Row count must match
#'       \code{nrow(X_train)}. Only applicable to \code{method = "dearseq"};
#'       ignored for other methods. Defaults to \code{NULL}, in which case
#'       the \code{covariates} argument is used for dearseq as well.}
#'     \item{\code{relative_gain_metric}}{Character string. Metric used to
#'       evaluate prediction quality. Only applicable to
#'       \code{method = "relative_gain"}. One of \code{"rmse"},
#'       \code{"srmse"}, \code{"r2"}, or \code{"spearman"}. Defaults to
#'       \code{"rmse"}. Ignored for other methods.}
#'     \item{\code{relative_gain_inner_folds}}{Positive integer. Number of
#'       inner CV folds. Only applicable to \code{method = "relative_gain"}.
#'       Defaults to \code{5}. Ignored for other methods.}
#'     \item{\code{relative_gain_seed}}{Integer. Random seed for inner fold
#'       assignment. Only applicable to \code{method = "relative_gain"}.
#'       Defaults to \code{12345}. Ignored for other methods.}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{selected_features}}{Character vector of selected column
#'       names, in decreasing order of score. For \code{"dearseq"} with
#'       \code{dearseq_level = "geneset"}, this is the union of member genes
#'       of the selected genesets (restricted to \code{colnames(X_train)}),
#'       not geneset names.}
#'     \item{\code{selection_method}}{Character string. The method used.}
#'     \item{\code{selection_scores}}{Named numeric vector of scores for
#'       \emph{all} features, in decreasing order. For
#'       \code{"relative_gain"}, scores are the gain values. For
#'       \code{"dearseq"} with \code{dearseq_level = "geneset"}, this is
#'       indexed by geneset name (not gene name). All \code{NA} when
#'       \code{top_n} exceeded (or equalled) the number of available
#'       features/genesets, since scores are not computed in that case.}
#'     \item{\code{top_n}}{Integer or \code{NULL}. The \code{top_n} value
#'       used after resolving precedence with \code{threshold} - capped to
#'       the number of available features/genesets if it exceeded that
#'       count.}
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
#'                 relative_gain_metric = "rmse",
#'                 relative_gain_inner_folds = 5))
#'
#' # Relative gain with covariates
#' covariates <- matrix(rnorm(40 * 2), nrow = 40,
#'                      dimnames = list(NULL, c("age", "sex")))
#' result <- run_selection(X, Y, covariates = covariates,
#'   params = list(method = "relative_gain", threshold = 0,
#'                 relative_gain_metric = "r2"))
#' }
#'
#' @export
# -----------------------------------------------------------------------------
run_selection <- function(X_train, Y_train = NULL, covariates = NULL,
                          treatment = NULL, individual_id = NULL,
                          timepoint = NULL, params) {

  # ---------------------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------------------
  .validate_X(X_train, allow_na = TRUE)
  .validate_selection_params(params, p = ncol(X_train))

  if (is.null(colnames(X_train)))
    stop("[predictomics] X_train must have column names for feature selection.",
         call. = FALSE)

  method    <- params$method
  top_n     <- params$top_n
  threshold <- params$threshold

  if (anyNA(X_train) && method %in% c("relative_gain", "rise", "dearseq"))
    stop(
      "[predictomics] selection_params$method = '", method, "' does not ",
      "support NA values in X_train. Please impute or remove missing ",
      "values first, or use a different method ('variance', 'pearson', ",
      "'spearman' support NA via na.rm/pairwise-complete observations).",
      call. = FALSE
    )

  # Supervised methods require Y_train
  if (method %in% c("pearson", "spearman", "relative_gain", "rise")) {
    if (is.null(Y_train))
      stop("[predictomics] Y_train must be provided for method = '",
           method, "'.", call. = FALSE)
    .validate_Y(Y_train)
    .validate_Y_X_compat(Y_train, X_train)
  }

  if (method == "rise") {

    rise_paired <- isTRUE(params$rise_paired %||% FALSE)

    if (rise_paired) {
      # Paired RISE is harmonised with dearseq's paired mode: the contrast
      # is timepoint (not treatment), and individual_id/timepoint pairing is
      # validated with the same shared helper.
      .validate_individual_timepoint_pairing(
        individual_id = individual_id,
        timepoint     = timepoint,
        n             = nrow(X_train),
        context       = "selection_params$method = 'rise' with rise_paired = TRUE"
      )
    } else {
      if (is.null(treatment))
        stop("[predictomics] treatment must be provided for method = 'rise'.",
             call. = FALSE)
      if (!is.numeric(treatment) || !all(treatment %in% c(0, 1)))
        stop("[predictomics] treatment must be a binary numeric vector (0/1) ",
             "for method = 'rise'.", call. = FALSE)
    }

    if (!requireNamespace("SurrogateRank", quietly = TRUE))
      stop("[predictomics] The SurrogateRank package is required for ",
           "method = 'rise'. Install it with: install.packages('SurrogateRank')",
           call. = FALSE)
  }

  if (method == "dearseq") {

    dearseq_mode <- params$dearseq_mode %||% "classic"

    if (dearseq_mode == "classic") {

      if (is.null(treatment))
        stop("[predictomics] treatment must be provided for method = ",
             "'dearseq' with dearseq_mode = 'classic'.", call. = FALSE)

      if (xor(is.null(individual_id), is.null(timepoint)))
        stop("[predictomics] For dearseq_mode = 'classic', 'individual_id' ",
             "and 'timepoint' must be supplied together, or not at all.",
             call. = FALSE)

      if (!is.null(individual_id))
        .validate_individual_timepoint_pairing(
          individual_id = individual_id,
          timepoint     = timepoint,
          n             = nrow(X_train),
          context       = paste0(
            "selection_params$method = 'dearseq' with dearseq_mode = ",
            "'classic' and individual_id/timepoint"
          )
        )

    } else {

      .validate_individual_timepoint_pairing(
        individual_id = individual_id,
        timepoint     = timepoint,
        n             = nrow(X_train),
        context       = "selection_params$dearseq_mode = 'paired'"
      )
    }

    dearseq_covariates <- params$dearseq_covariates
    if (!is.null(dearseq_covariates)) {
      if (!is.matrix(dearseq_covariates) && !is.data.frame(dearseq_covariates))
        stop("[predictomics] selection_params$dearseq_covariates must be a ",
             "numeric matrix or data frame, or NULL.", call. = FALSE)
      if (nrow(dearseq_covariates) != nrow(X_train))
        stop("[predictomics] selection_params$dearseq_covariates must have ",
             "the same number of rows as X_train (", nrow(X_train), ").",
             call. = FALSE)
      if (is.null(colnames(dearseq_covariates)))
        stop("[predictomics] selection_params$dearseq_covariates must have ",
             "column names.", call. = FALSE)
    }

    if (!requireNamespace("dearseq", quietly = TRUE))
      stop("[predictomics] The dearseq package is required for ",
           "method = 'dearseq'. Install it with: ",
           "BiocManager::install('dearseq')", call. = FALSE)
  }

  # top_n takes precedence - inform the user
  if (!is.null(top_n) && !is.null(threshold)) {
    message("[predictomics] Both top_n and threshold supplied to run_selection; ",
            "top_n takes precedence.")
    threshold <- NULL
  }

  # ---------------------------------------------------------------------------
  # 2. Compute scores (skipped entirely if top_n already covers every
  # available feature/geneset - selection scores would be discarded anyway
  # since all items are selected regardless of their value)
  # ---------------------------------------------------------------------------
  is_dearseq_geneset <- identical(method, "dearseq") &&
    identical(params$dearseq_level %||% "gene", "geneset") &&
    !is.null(params$genesets)
  effective_ids <- if (is_dearseq_geneset) names(params$genesets) else colnames(X_train)
  effective_p   <- length(effective_ids)
  unit_plural   <- if (is_dearseq_geneset) "genesets" else "features"

  top_n_exceeds_available <- !is.null(top_n) && top_n >= effective_p

  if (top_n_exceeds_available) {

    message(
      "[predictomics] selection_params$top_n (", top_n, ") exceeds (or ",
      "equals) the number of available ", unit_plural, " (", effective_p,
      "); selecting all ", effective_p, " ", unit_plural, " without ",
      "computing selection scores."
    )
    top_n  <- effective_p
    scores <- stats::setNames(rep(NA_real_, effective_p), effective_ids)

  } else {

    scores <- switch(method,

                     variance = {
                       apply(X_train, 2, var, na.rm = TRUE)
                     },

                     pearson = {
                       apply(X_train, 2, function(x)
                         abs(cor(x, Y_train, method = "pearson",
                                 use = "pairwise.complete.obs")))
                     },

                     spearman = {
                       apply(X_train, 2, function(x)
                         abs(cor(x, Y_train, method = "spearman",
                                 use = "pairwise.complete.obs")))
                     },

                     relative_gain = {
                       .compute_relative_gain(
                         X_train      = X_train,
                         Y_train      = Y_train,
                         covariates   = covariates,
                         metric       = params$relative_gain_metric       %||% "rmse",
                         inner_folds  = params$relative_gain_inner_folds  %||% 5L,
                         seed         = params$relative_gain_seed         %||% 12345L
                       )
                     },

                     rise = {
                       .compute_rise_scores(
                         X_train          = X_train,
                         Y_train          = Y_train,
                         treatment        = treatment,
                         individual_id    = individual_id,
                         timepoint        = timepoint,
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
                     },

                     dearseq = {
                       dearseq_covariates <- params$dearseq_covariates %||% covariates
                       if ((params$dearseq_level %||% "gene") == "geneset") {
                         .compute_dearseq_geneset_scores(
                           X_train         = X_train,
                           covariates      = dearseq_covariates,
                           treatment       = treatment,
                           individual_id   = individual_id,
                           timepoint       = timepoint,
                           dearseq_mode    = params$dearseq_mode          %||% "classic",
                           genesets        = params$genesets,
                           which_test      = params$dearseq_which_test    %||% "asymptotic",
                           preprocessed    = params$dearseq_preprocessed  %||% TRUE,
                           padjust_methods = params$dearseq_padjust_methods %||% "BH",
                           which_weights   = params$dearseq_which_weights %||% "loclin",
                           n_perm          = params$dearseq_n_perm        %||% 1000L,
                           bw              = params$dearseq_bw            %||% "nrd",
                           kernel          = params$dearseq_kernel        %||% "gaussian"
                         )
                       } else {
                         .compute_dearseq_scores(
                           X_train         = X_train,
                           covariates      = dearseq_covariates,
                           treatment       = treatment,
                           individual_id   = individual_id,
                           timepoint       = timepoint,
                           dearseq_mode    = params$dearseq_mode          %||% "classic",
                           which_test      = params$dearseq_which_test    %||% "asymptotic",
                           preprocessed    = params$dearseq_preprocessed  %||% TRUE,
                           padjust_methods = params$dearseq_padjust_methods %||% "BH",
                           which_weights   = params$dearseq_which_weights %||% "loclin",
                           n_perm          = params$dearseq_n_perm        %||% 1000L,
                           bw              = params$dearseq_bw            %||% "nrd",
                           kernel          = params$dearseq_kernel        %||% "gaussian"
                         )
                       }
                     }
    )

    if (method %in% c("rise", "dearseq")) {
      scores <- sort(scores, decreasing = FALSE)
    } else {
      scores <- sort(scores, decreasing = TRUE)
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Select features (or, for dearseq geneset-level, genesets - resolved to
  # gene names below)
  # ---------------------------------------------------------------------------
  unit <- if (is_dearseq_geneset) "geneset" else "feature"

  selected_labels <- if (!is.null(top_n)) {

    if (!top_n_exceeds_available && method == "relative_gain" &&
        any(scores[seq_len(top_n)] < 0)) {

      n_below_floor <- sum(scores[seq_len(top_n)] < 0)

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

    if (method %in% c("rise", "dearseq")) {

      # For RISE/dearseq, threshold is a maximum p-value (lower = better)
      sel <- names(scores)[scores <= threshold]

      if (length(sel) == 0L)
        stop("[predictomics] No ", unit, "s pass the '", method, "' p-value ",
             "threshold (", threshold, "). Consider raising the threshold or ",
             "using top_n instead.", call. = FALSE)

      sel

    } else {

      sel <- names(scores)[scores >= threshold]

      if (length(sel) == 0L)
        stop("[predictomics] No features pass the threshold (", threshold,
             ") for method = '", method, "'. ",
             "Consider lowering threshold or using top_n instead.",
             call. = FALSE)

      sel

    }
  }

  # For dearseq geneset-level, selected_labels are geneset names - resolve to
  # the union of their member genes (restricted to columns of X_train) so
  # that predict_cv()/run_fold() can subset X as usual.
  selected_features <- if (is_dearseq_geneset) {
    genes <- unique(unlist(params$genesets[selected_labels]))
    intersect(genes, colnames(X_train))
  } else {
    selected_labels
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
#' @param treatment Binary numeric vector (0/1). Treatment assignment. Used
#'   as the contrast when \code{paired = FALSE}; ignored when
#'   \code{paired = TRUE} (see \code{individual_id}/\code{timepoint}).
#' @param individual_id Vector identifying individuals, or \code{NULL}.
#'   Required when \code{paired = TRUE}: rows are sorted by
#'   \code{individual_id} (then \code{timepoint}) before splitting into
#'   pre/post groups, so that pre- and post-treatment rows for the same
#'   individual are matched regardless of input row order. Ignored when
#'   \code{paired = FALSE}.
#' @param timepoint Binary numeric vector (0/1), or \code{NULL}. Used as the
#'   contrast (0 = pre-treatment, 1 = post-treatment) when \code{paired =
#'   TRUE}, in place of \code{treatment}. Ignored when \code{paired = FALSE}.
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
                                 p_correction, n_cores, alternative, paired,
                                 individual_id = NULL, timepoint = NULL) {

  # ---------------------------------------------------------------------------
  # 1. Reshape inputs into RISE format
  # ---------------------------------------------------------------------------
  if (isTRUE(paired)) {

    # Sort by individual_id (then timepoint) so that pre-treatment
    # (timepoint == 0) and post-treatment (timepoint == 1) rows for the same
    # individual are matched by position, regardless of input row order.
    message(
      "[predictomics] RISE paired mode: rows are sorted by individual_id ",
      "(then timepoint) so that pre-treatment (timepoint == 0) and ",
      "post-treatment (timepoint == 1) rows for the same individual are ",
      "correctly matched, regardless of input row order."
    )

    ord       <- order(individual_id, timepoint)
    Y_ordered <- Y_train[ord]
    X_ordered <- X_train[ord, , drop = FALSE]
    tp_ordered <- timepoint[ord]

    idx1  <- which(tp_ordered == 1)
    idx0  <- which(tp_ordered == 0)

    yone  <- Y_ordered[idx1]
    yzero <- Y_ordered[idx0]
    sone  <- X_ordered[idx1, , drop = FALSE]
    szero <- X_ordered[idx0, , drop = FALSE]

  } else {

    idx1  <- which(treatment == 1)
    idx0  <- which(treatment == 0)

    yone  <- Y_train[idx1]
    yzero <- Y_train[idx0]
    sone  <- X_train[idx1, , drop = FALSE]
    szero <- X_train[idx0, , drop = FALSE]
  }

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
  # Internal scores for sorting (negated so higher = better)
  # Build interpretable score vector: adjusted p-value where < 1,
  # unadjusted p-value as tiebreaker for ceiling features
  scores_stored     <- p_adj
  at_ceiling        <- p_adj == 1
  if (any(at_ceiling))
    scores_stored[at_ceiling] <- p_unadj[at_ceiling]

  # scores_stored is on the raw p-value scale: smaller = more significant.
  # Returned unsorted; run_selection sorts ascending so best features come first.
  names(scores_stored) <- colnames(X_train)
  scores_stored
}


# -----------------------------------------------------------------------------
#' Compute dearseq adjusted p-value scores for feature selection
#'
#' @description
#' Calls \code{dearseq::dear_seq()} to test for differential expression
#' across a binary grouping variable, then returns adjusted p-values for all
#' features (columns) of \code{X_train}. Two modes are supported:
#' \describe{
#'   \item{\code{"classic"}}{Compares \code{treatment == 1} vs
#'     \code{treatment == 0}. If \code{individual_id}/\code{timepoint} are
#'     supplied, rows are first restricted to \code{timepoint == 1}.}
#'   \item{\code{"paired"}}{Compares \code{timepoint == 1} vs
#'     \code{timepoint == 0}, with \code{sample_group = individual_id}. If
#'     \code{treatment} is supplied, rows are first restricted to
#'     \code{treatment == 1}.}
#' }
#' Only rows are ever restricted before the comparison; all columns of
#' \code{X_train} are scored regardless.
#'
#' @param X_train Numeric matrix. Training predictor matrix (samples x
#'   features).
#' @param covariates Numeric matrix or data frame, or \code{NULL}. Passed to
#'   \code{dear_seq()} as a full design matrix (with intercept) via
#'   \code{\link{.build_dearseq_covariates}}.
#' @param treatment Binary numeric vector (0/1), factor, or \code{NULL}.
#' @param individual_id Vector identifying individuals, or \code{NULL}.
#' @param timepoint Binary numeric vector (0/1), or \code{NULL}.
#' @param dearseq_mode Character string. One of \code{"classic"} or
#'   \code{"paired"}.
#' @param which_test,preprocessed,padjust_methods,which_weights,n_perm,bw,kernel
#'   Arguments passed directly to \code{dearseq::dear_seq()}.
#'
#' @return Named numeric vector of adjusted p-values, one per feature, in
#'   \code{colnames(X_train)} order.
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_dearseq_scores <- function(X_train, covariates, treatment,
                                    individual_id, timepoint, dearseq_mode,
                                    which_test, preprocessed, padjust_methods,
                                    which_weights, n_perm, bw, kernel) {

  g <- .resolve_dearseq_groups(X_train, treatment, individual_id, timepoint,
                              dearseq_mode)

  X_sub   <- X_train[g$keep, , drop = FALSE]
  cov_sub <- if (!is.null(covariates)) covariates[g$keep, , drop = FALSE] else NULL

  exprmat <- t(X_sub)
  variables2test <- matrix(g$group_var, ncol = 1L,
                           dimnames = list(NULL, g$group_name))
  covariates_design <- .build_dearseq_covariates(cov_sub)

  res <- suppressMessages(
    dearseq::dear_seq(
      exprmat         = exprmat,
      covariates      = covariates_design,
      variables2test  = variables2test,
      sample_group    = g$sample_group,
      which_test      = which_test,
      preprocessed    = preprocessed,
      padjust_methods = padjust_methods,
      which_weights   = which_weights,
      n_perm          = n_perm,
      bw              = bw,
      kernel          = kernel
    )
  )

  adj_pval <- res[["pvals"]]$adjPval
  names(adj_pval) <- rownames(exprmat)

  adj_pval[colnames(X_train)]
}


# -----------------------------------------------------------------------------
#' Compute dearseq geneset-level adjusted p-value scores for feature selection
#'
#' @description
#' Calls \code{dearseq::dgsa_seq()} to test genesets (rather than individual
#' genes) for differential expression across a binary grouping variable, in
#' the same \code{"classic"}/\code{"paired"} row-restriction modes as
#' \code{\link{.compute_dearseq_scores}}. Returns adjusted p-values indexed
#' by geneset name (not gene name); \code{\link{run_selection}} resolves
#' selected genesets back to their member gene names.
#'
#' @param X_train Numeric matrix. Training predictor matrix (samples x
#'   features).
#' @param covariates Numeric matrix or data frame, or \code{NULL}. Passed to
#'   \code{dgsa_seq()} as a full design matrix (with intercept) via
#'   \code{\link{.build_dearseq_covariates}}.
#' @param treatment Binary numeric vector (0/1), factor, or \code{NULL}.
#' @param individual_id Vector identifying individuals, or \code{NULL}.
#' @param timepoint Binary numeric vector (0/1), or \code{NULL}.
#' @param dearseq_mode Character string. One of \code{"classic"} or
#'   \code{"paired"}.
#' @param genesets Named list of character vectors of feature names.
#' @param which_test,preprocessed,padjust_methods,which_weights,n_perm,bw,kernel
#'   Arguments passed directly to \code{dearseq::dgsa_seq()}.
#'
#' @return Named numeric vector of adjusted p-values, one per geneset, in
#'   \code{names(genesets)} order (restricted to genesets with at least one
#'   feature present in \code{X_train}).
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_dearseq_geneset_scores <- function(X_train, covariates, treatment,
                                            individual_id, timepoint,
                                            dearseq_mode, genesets, which_test,
                                            preprocessed, padjust_methods,
                                            which_weights, n_perm, bw, kernel) {

  g <- .resolve_dearseq_groups(X_train, treatment, individual_id, timepoint,
                              dearseq_mode)

  X_sub   <- X_train[g$keep, , drop = FALSE]
  cov_sub <- if (!is.null(covariates)) covariates[g$keep, , drop = FALSE] else NULL

  genesets_filtered <- .filter_genesets_to_available_features(
    X_sub, genesets, "dearseq (geneset-level)"
  )

  exprmat <- t(X_sub)
  variables2test <- matrix(g$group_var, ncol = 1L,
                           dimnames = list(NULL, g$group_name))
  covariates_design <- .build_dearseq_covariates(cov_sub)

  res <- suppressMessages(
    dearseq::dgsa_seq(
      exprmat         = exprmat,
      covariates      = covariates_design,
      variables2test  = variables2test,
      genesets        = genesets_filtered,
      sample_group    = g$sample_group,
      which_test      = which_test,
      preprocessed    = preprocessed,
      padjust_methods = padjust_methods,
      which_weights   = which_weights,
      n_perm          = n_perm,
      bw              = bw,
      kernel          = kernel
    )
  )

  adj_pval <- res[["pvals"]]$adjPval
  names(adj_pval) <- names(genesets_filtered)

  adj_pval
}


# -----------------------------------------------------------------------------
#' Resolve the grouping variable, sample_group, and row subset for dearseq
#'
#' @description
#' Shared row-restriction logic for \code{\link{.compute_dearseq_scores}} and
#' \code{\link{.compute_dearseq_geneset_scores}}. For \code{"classic"} mode,
#' restricts to \code{timepoint == 1} rows when \code{individual_id} is
#' supplied and returns \code{treatment} as the grouping variable. For
#' \code{"paired"} mode, restricts to \code{treatment == 1} rows when
#' \code{treatment} is supplied and returns \code{timepoint} as the grouping
#' variable with \code{individual_id} as \code{sample_group}.
#'
#' @param X_train Numeric matrix. Training predictor matrix.
#' @param treatment Binary numeric vector (0/1), factor, or \code{NULL}.
#' @param individual_id Vector identifying individuals, or \code{NULL}.
#' @param timepoint Binary numeric vector (0/1), or \code{NULL}.
#' @param dearseq_mode Character string. One of \code{"classic"} or
#'   \code{"paired"}.
#'
#' @return A named list with \code{keep} (integer row indices), \code{group_var}
#'   (numeric grouping vector), \code{group_name} (character), and
#'   \code{sample_group} (vector or \code{NULL}).
#' @keywords internal
# -----------------------------------------------------------------------------
.resolve_dearseq_groups <- function(X_train, treatment, individual_id,
                                    timepoint, dearseq_mode) {

  treatment_bin <- if (!is.null(treatment)) .coerce_treatment_binary(treatment)
  else NULL

  if (dearseq_mode == "classic") {

    keep <- if (!is.null(individual_id))
      which(timepoint == 1)
    else
      seq_len(nrow(X_train))

    list(
      keep         = keep,
      group_var    = treatment_bin[keep],
      group_name   = "treatment",
      sample_group = NULL
    )

  } else {

    keep <- if (!is.null(treatment_bin))
      which(treatment_bin == 1)
    else
      seq_len(nrow(X_train))

    list(
      keep         = keep,
      group_var    = timepoint[keep],
      group_name   = "timepoint",
      sample_group = individual_id[keep]
    )
  }
}
