# =============================================================================
# predict_cv.R
# Master cross-validation prediction function for the predictomics package.
# =============================================================================


# -----------------------------------------------------------------------------
#' Cross-validated prediction from gene expression data
#'
#' @description
#' Master function for end-to-end prediction from gene expression data in a
#' cross-validation (CV) framework. Optionally applies feature engineering,
#' feature selection, and a choice of predictive model within each CV fold,
#' returning observed and predicted values alongside fold assignments and
#' pipeline metadata.
#'
#' @details
#' The pipeline is applied sequentially in the order:
#' feature engineering -> feature selection -> model fitting -> prediction.
#' By default, all steps are performed **within** each CV fold, meaning that
#' engineering and selection are re-fitted on the training partition and then
#' applied to the held-out test partition. This is the statistically correct
#' approach and avoids data leakage.
#'
#' If \code{outside_cv = TRUE}, engineering and selection steps are instead
#' applied to the full dataset before the CV loop begins.
#' **This will produce optimistically biased performance estimates and is
#' encouraged to be used for exploratory analysis only.**
#'
#' If \code{engineering_params$gene_level_fc = TRUE} is specified,
#' \code{individual_id} and \code{timepoint} must also be supplied. Every
#' individual must have exactly two observations, one at
#' \code{timepoint = 0} (pre-treatment) and one at \code{timepoint = 1}
#' (post-treatment); otherwise an informative error is raised. Before the CV
#' loop, any other engineering, or selection is performed, each individual's
#' paired rows are replaced by a single row equal to the post-treatment minus
#' pre-treatment value for every feature (and the response \code{Y} is taken
#' from the post-treatment observation). \code{individual_id} and
#' \code{timepoint} are then discarded and the remainder of the pipeline
#' operates on the fold-change-transformed dataset as normal. This mode is
#' not compatible with \code{selection_params$rise_paired = TRUE}.
#'
#' If \code{selection_params$method = "dearseq"} is specified, a
#' differential expression filter (\code{dearseq::dear_seq()}, or
#' \code{dearseq::dgsa_seq()} when
#' \code{selection_params$dearseq_level = "geneset"}) is applied
#' \strong{once, on the full, untransformed dataset, before any engineering
#' (including \code{gene_level_fc}) or other selection, and regardless of
#' \code{outside_cv}} - differential expression analysis is not meaningful on
#' transformed data, so this is a deliberate exception to the fold-safe
#' default described above. The resulting gene list is used to filter
#' \code{X} before the rest of the pipeline runs as normal; \code{dearseq}
#' fully replaces the selection step for that call (no further per-fold or
#' outside-cv selection is performed). See \code{\link{run_selection}} for
#' the \code{"classic"}/\code{"paired"} modes, the \code{"gene"}/
#' \code{"geneset"} levels, and hyperparameters.
#' \code{selection_params$dearseq_level = "gene"} (the default) cannot be
#' combined with \code{engineering_params$genesets}, since gene-level DEA
#' filtering typically leaves most genesets with no surviving member genes;
#' use \code{dearseq_level = "geneset"} instead if geneset-level engineering
#' is required downstream. With \code{dearseq_level = "geneset"}, any
#' geneset in \code{engineering_params$genesets} left with no surviving
#' genes after the dearseq filter is automatically dropped before
#' engineering runs (with a message, if \code{verbose = TRUE}). If
#' \code{selection_params$dearseq_mode = "paired"} is also specified, this
#' upfront differential expression filter is a \strong{filtration step only}:
#' it uses both pre- and post-treatment rows to score genes, but once genes
#' are selected, the dataset is immediately restricted to post-treatment
#' (\code{timepoint == 1}) rows for every subsequent stage (engineering,
#' selection, model fitting, prediction); pre-treatment rows are discarded
#' and play no further role. This exactly mirrors how
#' \code{selection_params$rise_paired = TRUE} handles its own paired
#' filtration step below, so the two paired approaches behave identically
#' with respect to which rows continue on to modelling. This mode is not
#' compatible with \code{engineering_params$gene_level_fc = TRUE}.
#'
#' If \code{selection_params$method = "rise"} with
#' \code{selection_params$rise_paired = TRUE} is specified, \code{individual_id}
#' and \code{timepoint} must also be supplied, with the same pairing
#' requirement as \code{gene_level_fc}/\code{dearseq_mode = "paired"} (every
#' individual has exactly one \code{timepoint = 0} and one
#' \code{timepoint = 1} observation), validated by the same shared check.
#' \code{timepoint} is used as the pre/post contrast for RISE screening, in
#' place of \code{treatment} (which, if supplied, retains its usual role as
#' an optional model predictor via \code{treatment_predictor}, decoupled from
#' the RISE pairing). Rows are sorted by \code{individual_id} (then
#' \code{timepoint}) before fitting, so pre- and post-treatment rows for the
#' same individual are matched correctly regardless of input row order. As
#' with \code{dearseq_mode = "paired"} above, RISE screening (the filtration
#' step) uses both pre- and post-treatment rows, but modelling is then
#' performed on post-treatment (\code{timepoint == 1}) rows only;
#' pre-treatment rows are discarded after screening. This mode is not
#' compatible with \code{engineering_params$gene_level_fc = TRUE}. Unlike
#' \code{selection_params$method = "dearseq"}, RISE has no geneset-level
#' counterpart: paired RISE screening always operates on the raw, gene-level
#' predictor matrix, so it is also not compatible with
#' \code{engineering_params$genesets} (geneset-level engineering).
#'
#' If \code{selection_params$method = "variance"} is combined with
#' \code{engineering_params$col_transform = "z"}, a warning is issued and
#' variance scores are computed on a separately re-engineered, pre-z-score
#' version of the training data instead (identical geneset aggregation, if
#' \code{engineering_params$genesets} is set, just without the z-score step).
#' This is necessary because z-scoring fixes every feature's own variance to
#' exactly ~1 on the training fold, before any aggregation - so filtering by
#' variance on the fully engineered matrix would select on essentially
#' floating-point noise (with no aggregation), or on a matrix confounded by
#' geneset size and internal gene-gene correlation rather than genuine
#' cross-sample variability (with aggregation). Only the variance
#' \emph{selection scores} are computed differently; the training/test
#' matrices actually used for model fitting still come from the normal,
#' fully engineered (z-scored) pipeline, subset to whichever features
#' variance selection chose. No other selection method is affected by this,
#' since \code{"pearson"}/\code{"spearman"} correlation with \code{Y} is
#' invariant to each feature's own linear rescaling.
#'
#' If \code{treatment} is supplied and \code{treatment_predictor = TRUE},
#' treatment is appended to the predictor matrix after engineering and
#' selection, immediately before model fitting. Treatment is never passed
#' through engineering or selection. If \code{treatment} is a factor with
#' k levels, k-1 dummy columns are created using the first level as reference
#' (via \code{model.matrix}). If \code{treatment} is binary numeric, it is
#' appended as a single column named \code{"treatment"}.
#'
#' If \code{covariates} is supplied, it is one-hot encoded once before the CV
#' loop and appended to the predictor matrix after engineering, selection, and
#' treatment, immediately before model fitting. Covariates are never passed
#' through engineering or selection. Covariate column names must not clash with
#' feature names in \code{X} or treatment column names.
#'
#' If \code{X = NULL}, a **baseline model** is fit instead of the usual
#' gene-expression pipeline: engineering and selection are skipped entirely
#' (both \code{engineering_params} and \code{selection_params} must be
#' \code{NULL}), and each fold's model is fit using only \code{covariates}
#' and/or \code{treatment} (if \code{treatment_predictor = TRUE}) as
#' predictors. If neither is available either, the "model" is simply the
#' training-fold mean of \code{Y}, predicted for every held-out sample - no
#' call to \code{\link{run_model}} is made in that case. An informative
#' message is printed (when \code{verbose = TRUE}) indicating that a
#' baseline model is being used, and whether it has any predictors at all.
#'
#' **Parallelisation**: the outer CV loop is parallelised via
#' \code{future.apply::future_lapply}. To enable parallel execution, set a
#' \code{future} plan before calling \code{predict_cv}:
#' \preformatted{
#'   future::plan(future::multisession, workers = 4)
#'   result <- predict_cv(Y, X, ...)
#'   future::plan(future::sequential)  # reset after use
#' }
#' If no plan is set, execution is sequential (the default behaviour).
#' Reproducibility is ensured via \code{future.seed} regardless of whether
#' parallel or sequential execution is used. Per-fold progress messages are
#' suppressed when a parallel plan is active.
#'
#' @param Y Numeric vector of length n. The response variable to be predicted.
#' @param X Numeric matrix of dimensions n (samples) x p (features).
#'   The predictor matrix. Column names should be feature identifiers. Pass
#'   \code{NULL} (default) to fit a \strong{baseline model} instead: only
#'   \code{covariates} (and \code{treatment}, if \code{treatment_predictor =
#'   TRUE}) are used as predictors, or, if none of those are supplied either,
#'   the model simply predicts the training-fold mean of \code{Y}. When
#'   \code{X = NULL}, \code{engineering_params} and \code{selection_params}
#'   must also be \code{NULL}, since they operate on gene-level features.
#'   \code{X} may contain \code{NA} values; support depends on the chosen
#'   \code{engineering_params$col_transform}/\code{agg_method},
#'   \code{selection_params$method}, and \code{model_params$impute} - see
#'   \code{\link{run_engineering}}, \code{\link{run_selection}}, and
#'   \code{\link{run_model}} for details on which options tolerate \code{NA}.
#' @param cv_type Character string. Type of cross-validation. One of
#'   \code{"kfold"} (K-fold CV) or \code{"loo"} (leave-one-out CV).
#'   Defaults to \code{"kfold"}.
#' @param folds Positive integer. Number of folds for k-fold CV. Must satisfy
#'   \code{2 <= folds <= n}. Ignored when \code{cv_type = "loo"}.
#'   Defaults to \code{10}.
#' @param seed Integer. Random seed for reproducible fold assignment and
#'   parallel RNG streams. Defaults to \code{12345}.
#' @param engineering_params A named list specifying feature engineering steps
#'   to apply. See \code{\link{run_engineering}} for supported options. Pass
#'   \code{NULL} (default) to skip engineering.
#' @param selection_params A named list specifying feature selection steps to
#'   apply. See \code{\link{run_selection}} for supported options. Pass
#'   \code{NULL} (default) to skip selection.
#' @param model_params A named list specifying the model and hyperparameter
#'   options. \code{model_params$method} must be one of \code{"lm"}
#'   (ordinary least squares), \code{"glmnet"} (elastic net), \code{"ridge"},
#'   \code{"lasso"}, \code{"ranger"} (random forest), or \code{"svr"}
#'   (support vector regression, linear kernel). See \code{\link{run_model}}
#'   for supported options. Defaults to \code{list(method = "lm")}
#'   (ordinary least squares).
#' @param outside_cv Logical. If \code{TRUE}, engineering and selection steps
#'   are applied to the full dataset \strong{before} the CV loop. This
#'   introduces data leakage and may produce biased performance estimates.
#'   Defaults to \code{FALSE}.
#' @param treatment A factor or binary numeric vector of length n encoding
#'   treatment group membership. If a factor, levels are used as group labels.
#'   If numeric, must contain only 0 and 1. Pass \code{NULL} (default) if no
#'   treatment variable is available.
#' @param treatment_predictor Logical. If \code{TRUE}, treatment is appended
#'   to the predictor matrix (after engineering and selection) as a covariate.
#'   Ignored when \code{treatment = NULL}. Defaults to \code{FALSE}.
#' @param covariates A numeric matrix or data frame of dimensions n (samples)
#'   x q (covariates). Covariates are protected predictors always included in
#'   model fitting. Factor columns are one-hot encoded via
#'   \code{model.matrix}. Pass \code{NULL} (default) for no covariates.
#' @param individual_id A vector of length n identifying which individual each
#'   sample (row of \code{X}) belongs to. Required when
#'   \code{engineering_params$gene_level_fc = TRUE} is used, in which case
#'   every individual must have exactly two observations: one at
#'   \code{timepoint = 0} (pre-treatment) and one at \code{timepoint = 1}
#'   (post-treatment). Also required when \code{selection_params$dearseq_mode
#'   = "paired"} or \code{selection_params$rise_paired = TRUE} is used (same
#'   pairing requirement), and optionally used to restrict
#'   \code{selection_params$dearseq_mode = "classic"} to
#'   \code{timepoint == 1} rows. Pass \code{NULL} (default) otherwise.
#' @param timepoint A binary numeric vector of length n (0 = pre-treatment,
#'   1 = post-treatment), paired with \code{individual_id}. Required when
#'   \code{engineering_params$gene_level_fc = TRUE},
#'   \code{selection_params$dearseq_mode = "paired"}, or
#'   \code{selection_params$rise_paired = TRUE} is used - for the latter,
#'   \code{timepoint} is the pre/post contrast (in place of \code{treatment})
#'   and rows are sorted by \code{individual_id} (then \code{timepoint})
#'   before fitting, so input row order does not matter. Pass \code{NULL}
#'   (default) otherwise.
#' @param verbose Logical. If \code{TRUE}, prints progress messages throughout.
#'   Defaults to \code{TRUE}.
#'
#' @return An object of class \code{"predictomics"}, which is a named list
#'   containing:
#'   \describe{
#'     \item{\code{observed}}{Numeric vector of observed \code{Y} values, in
#'       original sample order.}
#'     \item{\code{predicted}}{Numeric vector of cross-validated predicted
#'       values, in original sample order.}
#'     \item{\code{fold_ids}}{Integer vector of fold assignments (1 to
#'       \code{folds}), in original sample order.}
#'     \item{\code{treatment}}{The \code{treatment} argument as supplied, or
#'       \code{NULL}. Restricted to post-treatment (\code{timepoint == 1})
#'       rows, matching \code{observed}/\code{predicted}, when
#'       \code{selection_params$rise_paired = TRUE} or
#'       \code{selection_params$dearseq_mode = "paired"} is used.}
#'     \item{\code{treatment_predictor}}{Logical. Whether treatment was used
#'       as a predictor.}
#'     \item{\code{covariates}}{The \code{covariates} argument as supplied, or
#'       \code{NULL}.}
#'     \item{\code{engineering_params}}{The \code{engineering_params} argument
#'       as supplied.}
#'     \item{\code{selection_params}}{The \code{selection_params} argument as
#'       supplied.}
#'     \item{\code{model_params}}{The \code{model_params} argument as
#'       supplied.}
#'     \item{\code{outside_cv}}{Logical. Whether outside-CV mode was used.}
#'     \item{\code{cv_type}}{Character string. The CV type used.}
#'     \item{\code{n_folds}}{Integer. Number of folds used.}
#'     \item{\code{n_samples}}{Integer. Number of samples. For
#'       \code{selection_params$rise_paired = TRUE} or
#'       \code{selection_params$dearseq_mode = "paired"}, this is the
#'       pre-restriction (both-arms) count, same as \code{n_samples_total};
#'       see \code{n_samples_modelled} for the post-treatment-only count
#'       actually used for modelling.}
#'     \item{\code{n_samples_modelled}}{Integer. Number of samples actually
#'       used for modelling (i.e. \code{length(observed)}). Equal to
#'       \code{n_samples} unless \code{selection_params$rise_paired = TRUE}
#'       or \code{selection_params$dearseq_mode = "paired"} restricted
#'       modelling to post-treatment rows only.}
#'     \item{\code{n_samples_total}}{Integer. Same as \code{n_samples}.}
#'     \item{\code{paired_rise}}{Logical. Whether
#'       \code{selection_params$rise_paired = TRUE} was used.}
#'     \item{\code{dearseq_paired}}{Logical. Whether
#'       \code{selection_params$method = "dearseq"} with
#'       \code{dearseq_mode = "paired"} was used.}
#'     \item{\code{n_features_input}}{Integer. Number of features in the input
#'       \code{X}. \code{0} if \code{X = NULL} (baseline model).}
#'     \item{\code{baseline_model}}{Logical. Whether a baseline model
#'       (\code{X = NULL}) was fit.}
#'     \item{\code{dearseq_selection}}{A list containing
#'       \code{selected_features} (gene names for
#'       \code{dearseq_level = "gene"}; \strong{geneset} names for
#'       \code{dearseq_level = "geneset"} - \code{X} itself is still
#'       filtered by the resolved member genes of those genesets),
#'       \code{selection_scores} (adjusted p-values - per gene for
#'       \code{dearseq_level = "gene"}, per geneset for
#'       \code{dearseq_level = "geneset"}), \code{dearseq_mode},
#'       \code{dearseq_level}, and \code{n_selected} (number of genes or
#'       genesets selected, matching \code{dearseq_level}) from the upfront
#'       dearseq filter. \code{NULL} if
#'       \code{selection_params$method != "dearseq"}.}
#'     \item{\code{fold_selection_diagnostics}}{A list of length
#'       \code{n_folds}, each element containing \code{selected_features},
#'       \code{selection_scores}, and \code{n_selected} for that fold.
#'       \code{NULL} if no explicit selection was performed.}
#'     \item{\code{outside_cv_selection}}{A list containing
#'       \code{selected_features}, \code{selection_scores}, and
#'       \code{n_selected} from outside-CV selection. \code{NULL} if
#'       \code{outside_cv = FALSE} or no selection was performed.}
#'     \item{\code{fold_embedded_selection_diagnostics}}{A list of length
#'       \code{n_folds}, each element containing \code{selected_features}
#'       and \code{n_selected} from embedded selection (lasso/glmnet).
#'       \code{NULL} if the model does not perform embedded selection.}
#'     \item{\code{fold_feature_importance}}{A list of length \code{n_folds},
#'       each element containing \code{scores} (named numeric vector) and
#'       \code{type} (\code{"coefficient"} or \code{"impurity"}) from that
#'       fold's fitted model. \code{NULL} unless
#'       \code{model_params$compute_importance = TRUE}. See
#'       \code{\link{run_model}} and \code{\link{plot_feature_importance}}.}
#'     \item{\code{call}}{The matched call.}
#'   }
#'
#' @seealso
#'   \code{\link{make_folds}} for fold generation,
#'   \code{\link{run_engineering}} for feature engineering options,
#'   \code{\link{run_selection}} for feature selection options,
#'   \code{\link{run_model}} for model options,
#'   \code{\link{print.predictomics}},
#'   \code{\link{plot.predictomics}},
#'   \code{\link{metrics.predictomics}} for result inspection.
#'
#' @examples
#' \dontrun{
#' # Simulate data: 50 samples, 200 genes
#' set.seed(1)
#' n <- 50; p <- 200
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' colnames(X) <- paste0("gene", seq_len(p))
#' rownames(X) <- paste0("sample", seq_len(n))
#' Y <- X[, 1] * 2 + rnorm(n)
#'
#' # Basic usage: 10-fold CV with default OLS model
#' result <- predict_cv(Y = Y, X = X)
#'
#' # With parallelisation
#' future::plan(future::multisession, workers = 4)
#' result <- predict_cv(Y = Y, X = X)
#' future::plan(future::sequential)
#'
#' # With binary treatment as predictor
#' treatment <- sample(c(0L, 1L), n, replace = TRUE)
#' result <- predict_cv(Y = Y, X = X, treatment = treatment,
#'                      treatment_predictor = TRUE)
#'
#' # With covariates
#' covariates <- data.frame(age = rnorm(n), sex = factor(sample(c("M","F"), n,
#'                          replace = TRUE)))
#' result <- predict_cv(Y = Y, X = X, covariates = covariates)
#'
#' # Baseline model: X = NULL, predicting from covariates alone
#' result <- predict_cv(Y = Y, X = NULL, covariates = covariates)
#'
#' # Null baseline model: X = NULL and no covariates - predicts the
#' # training-fold mean of Y
#' result <- predict_cv(Y = Y, X = NULL)
#' }
#'
#' @export
# -----------------------------------------------------------------------------
predict_cv <- function(Y,
                       X                   = NULL,
                       cv_type             = "kfold",
                       folds               = 10L,
                       seed                = 12345L,
                       engineering_params  = NULL,
                       selection_params    = NULL,
                       model_params        = list(method = "lm"),
                       outside_cv          = FALSE,
                       treatment           = NULL,
                       treatment_predictor = FALSE,
                       covariates          = NULL,
                       individual_id       = NULL,
                       timepoint           = NULL,
                       verbose             = TRUE) {

  cl <- match.call()

  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  n <- length(Y)
  is_baseline_model <- is.null(X)
  p <- if (is_baseline_model) 0L else ncol(X)
  folds <- if (cv_type == "loo") n else folds

  .validate_Y(Y)
  if (!is_baseline_model) {
    .validate_X(X, allow_na = TRUE)
    .validate_Y_X_compat(Y, X)
  }
  .validate_scalar_args(cv_type = cv_type, folds = folds, n = n,
                        seed = seed, outside_cv = outside_cv,
                        verbose = verbose,
                        treatment_predictor = treatment_predictor)
  .validate_treatment(treatment, Y)
  .validate_covariates(covariates, Y)
  .validate_params_list(params = engineering_params, arg_name = "engineering_params")
  .validate_params_list(params = selection_params,   arg_name = "selection_params")
  .validate_params_list(params = model_params,       arg_name = "model_params")

  # ---------------------------------------------------------------------------
  # 1a-baseline. X = NULL: fit a baseline model using only covariates/treatment
  # (or, if neither is available, the training-fold mean of Y). Gene-level
  # engineering/selection require X and are therefore disallowed here.
  # ---------------------------------------------------------------------------
  if (is_baseline_model) {

    if (!is.null(engineering_params) || !is.null(selection_params))
      stop(
        "[predictomics] X = NULL (baseline model) is not compatible with ",
        "engineering_params or selection_params, since these operate on ",
        "gene-level features. Set both to NULL (the default) to fit a ",
        "baseline model.",
        call. = FALSE
      )

    X <- matrix(nrow = n, ncol = 0)

    baseline_has_predictors <- !is.null(covariates) ||
      (!is.null(treatment) && treatment_predictor)

    if (verbose) {
      if (baseline_has_predictors) {
        message(
          "[predictomics] X = NULL: fitting a baseline model using only ",
          "covariates and/or treatment as predictors (no gene-level ",
          "features)."
        )
      } else {
        message(
          "[predictomics] X = NULL, with no covariates and no treatment ",
          "predictor: fitting a null baseline model that simply predicts ",
          "the training-fold mean of Y."
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 1a. Detect gene-level fold-change (gene_level_fc) engineering
  # ---------------------------------------------------------------------------
  is_gene_level_fc <- !is.null(engineering_params) &&
    isTRUE(engineering_params$gene_level_fc)

  # treatment_pipeline/covariates_pipeline/selection_params_pipeline/
  # engineering_params_pipeline are used for all downstream computation;
  # treatment/covariates/selection_params/engineering_params are left
  # untouched so that the returned result object reflects the arguments
  # exactly as supplied by the user.
  treatment_pipeline        <- treatment
  covariates_pipeline       <- covariates
  selection_params_pipeline <- selection_params
  engineering_params_pipeline <- engineering_params

  # ---------------------------------------------------------------------------
  # 1b. Detect paired RISE mode
  # ---------------------------------------------------------------------------
  is_paired_rise <- !is.null(selection_params) &&
    isTRUE(selection_params$method == "rise") &&
    isTRUE(selection_params$rise_paired)

  if (is_gene_level_fc && is_paired_rise) {
    stop(
      "[predictomics] engineering_params$gene_level_fc = TRUE is not ",
      "compatible with selection_params$rise_paired = TRUE. Set ",
      "rise_paired = FALSE (or omit it) to use gene_level_fc.",
      call. = FALSE
    )
  }

  # Paired RISE always screens the raw, un-engineered gene-level matrix (see
  # run_fold()'s is_paired_rise branch, which scores X_full rather than the
  # already-engineered X_train), so its selected feature names are gene
  # names. If engineering_params$genesets is set, the (geneset-named)
  # engineered matrix would then be subset by those gene names, which do not
  # match any column - this fails with an uninformative "subscript out of
  # bounds" error deep inside run_fold() rather than here, so it is rejected
  # explicitly upfront instead. Unlike gene_level_fc/dearseq_mode = "paired",
  # RISE has no geneset-level counterpart to fall back to.
  if (is_paired_rise && !is.null(engineering_params) &&
      !is.null(engineering_params$genesets)) {
    stop(
      "[predictomics] selection_params$rise_paired = TRUE is not ",
      "compatible with engineering_params$genesets (geneset-level ",
      "engineering). RISE paired screening always operates on the raw, ",
      "gene-level predictor matrix, so its selected genes cannot be ",
      "resolved against geneset-aggregated features. Remove ",
      "engineering_params$genesets, or use a different selection method.",
      call. = FALSE
    )
  }

  # is_dearseq / is_dearseq_paired are determined here (ahead of the dearseq
  # filtering block below) purely so the gene_level_fc incompatibility check
  # can live alongside the RISE one above.
  is_dearseq <- !is.null(selection_params) &&
    isTRUE(selection_params$method == "dearseq")
  is_dearseq_paired <- is_dearseq &&
    identical(selection_params$dearseq_mode %||% "classic", "paired")

  if (is_gene_level_fc && is_dearseq_paired) {
    stop(
      "[predictomics] engineering_params$gene_level_fc = TRUE is not ",
      "compatible with selection_params$method = 'dearseq' with ",
      "dearseq_mode = 'paired'. Set dearseq_mode = 'classic' (or use a ",
      "different selection method) to use gene_level_fc.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 1b2. Apply dearseq upfront gene filtering (before CV, before any
  # engineering - including gene_level_fc - or other selection). Runs on the
  # raw, untransformed data regardless of outside_cv, since differential
  # expression analysis is not meaningful on transformed data.
  # ---------------------------------------------------------------------------
  dearseq_selection <- NULL

  if (is_dearseq) {

    dearseq_level <- selection_params$dearseq_level %||% "gene"

    if (dearseq_level == "gene" && !is.null(engineering_params) &&
        !is.null(engineering_params$genesets)) {
      stop(
        "[predictomics] selection_params$method = 'dearseq' with ",
        "dearseq_level = 'gene' is not compatible with ",
        "engineering_params$genesets (geneset-level engineering). ",
        "Gene-level DEA filtering typically leaves most genesets with no ",
        "surviving member genes, causing downstream errors. Either set ",
        "selection_params$dearseq_level = 'geneset' to filter by ",
        "differentially expressed genesets instead, or remove ",
        "engineering_params$genesets.",
        call. = FALSE
      )
    }

    if (verbose) {
      message(
        "[predictomics] Applying dearseq differential expression filtering ",
        "once on the full, untransformed dataset (not per fold, and ",
        "regardless of outside_cv), since DEA is not meaningful on ",
        "transformed data. The resulting gene list is then used to filter ",
        "the dataset before any further engineering or selection."
      )
    }

    dea_fit <- run_selection(
      X_train       = X,
      Y_train       = NULL,
      covariates    = covariates,
      treatment     = treatment,
      individual_id = individual_id,
      timepoint     = timepoint,
      params        = selection_params
    )

    X <- X[, dea_fit$selected_features, drop = FALSE]
    p <- ncol(X)

    # For dearseq_level = "geneset", report selection at the geneset level
    # (selected genesets, not their resolved member genes) in the returned
    # diagnostics, even though X is filtered by the resolved gene names above.
    dearseq_selected <- if (dearseq_level == "geneset") {
      if (!is.null(dea_fit$top_n)) {
        names(dea_fit$selection_scores)[seq_len(dea_fit$top_n)]
      } else {
        names(dea_fit$selection_scores)[dea_fit$selection_scores <= dea_fit$threshold]
      }
    } else {
      dea_fit$selected_features
    }

    dearseq_selection <- list(
      selected_features = dearseq_selected,
      selection_scores  = dea_fit$selection_scores,
      dearseq_mode      = selection_params$dearseq_mode %||% "classic",
      dearseq_level     = dearseq_level,
      n_selected        = length(dearseq_selected)
    )

    # dearseq fully replaces the selection step for this call - suppress the
    # normal per-fold/outside-cv selection logic downstream.
    selection_params_pipeline <- NULL

    # If engineering_params$genesets is specified (only possible here when
    # dearseq_level = "geneset", since dearseq_level = "gene" is rejected
    # above), some genesets may have no genes left in the dearseq-filtered X.
    # Drop them so downstream engineering doesn't error on empty overlaps.
    if (!is.null(engineering_params) && !is.null(engineering_params$genesets)) {

      genesets_filtered <- lapply(engineering_params$genesets, intersect,
                                  y = colnames(X))
      keep    <- vapply(genesets_filtered, length, integer(1)) > 0L
      removed <- names(engineering_params$genesets)[!keep]

      if (length(removed) > 0L) {
        if (verbose) {
          message(
            "[predictomics] Removed ", length(removed), " geneset(s) from ",
            "engineering_params$genesets with no surviving genes after the ",
            "dearseq filter: ", paste(removed, collapse = ", "), "."
          )
        }
        engineering_params_pipeline$genesets <-
          engineering_params$genesets[keep]
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 1b3. dearseq_mode = "paired": after the filtration step above (which used
  # both pre- and post-treatment rows), restrict the dataset to post-treatment
  # (timepoint == 1) rows for the modelling stage, discarding pre-treatment
  # rows. This mirrors selection_params$rise_paired = TRUE's row handling
  # exactly, so that both paired filtration approaches leave the same rows
  # available for engineering/selection/modelling downstream.
  # ---------------------------------------------------------------------------
  n_samples_total_paired <- NULL

  if (is_dearseq_paired) {

    .validate_individual_timepoint_pairing(
      individual_id = individual_id,
      timepoint     = timepoint,
      n             = n,
      context       = "selection_params$dearseq_mode = 'paired'"
    )

    if (verbose) {
      message(
        "[predictomics] dearseq_mode = 'paired': the differential ",
        "expression filtering above used both pre- and post-treatment ",
        "rows. Predictive modelling will be performed on post-treatment ",
        "samples only (timepoint == 1); pre-treatment samples are ",
        "discarded from the modelling stage."
      )
    }

    keep <- which(timepoint == 1)

    n_samples_total_paired <- n

    X                   <- X[keep, , drop = FALSE]
    Y                   <- Y[keep]
    treatment_pipeline  <- if (!is.null(treatment_pipeline))
      treatment_pipeline[keep] else NULL
    covariates_pipeline <- if (!is.null(covariates_pipeline))
      covariates_pipeline[keep, , drop = FALSE] else NULL
    individual_id       <- individual_id[keep]
    timepoint           <- timepoint[keep]
    n                   <- length(Y)

    if (folds > n)
      stop(
        "[predictomics] folds (", folds, ") exceeds the number of ",
        "post-treatment samples (", n, ") after dearseq paired filtering. ",
        "Reduce 'folds' or use cv_type = 'loo'.",
        call. = FALSE
      )
    folds <- if (cv_type == "loo") n else folds
  }

  # ---------------------------------------------------------------------------
  # 1c. Apply gene-level fold-change transformation (before CV, before any
  # other engineering/selection)
  # ---------------------------------------------------------------------------
  if (is_gene_level_fc) {

    .validate_individual_timepoint_pairing(individual_id = individual_id,
                                           timepoint     = timepoint,
                                           n             = length(Y))

    fc_fit <- .compute_gene_level_fc(X             = X,
                                     Y             = Y,
                                     individual_id = individual_id,
                                     timepoint     = timepoint,
                                     treatment     = treatment,
                                     covariates    = covariates)

    if (verbose) {
      message(
        "[predictomics] Applied gene_level_fc engineering: collapsed ", n,
        " paired observations into ", nrow(fc_fit$X), " individual-level ",
        "fold-change rows (post-treatment - pre-treatment)."
      )
    }

    X                   <- fc_fit$X
    Y                   <- fc_fit$Y
    treatment_pipeline  <- fc_fit$treatment
    covariates_pipeline <- fc_fit$covariates
    n                   <- length(Y)
    p                   <- ncol(X)

    if (folds > n)
      stop(
        "[predictomics] folds (", folds, ") exceeds the number of ",
        "individuals (", n, ") after the gene_level_fc transformation. ",
        "Reduce 'folds' or use cv_type = 'loo'.",
        call. = FALSE
      )
    folds <- if (cv_type == "loo") n else folds
  }

  if (is_paired_rise) {

    # Harmonised with dearseq's paired mode: the contrast is timepoint (not
    # treatment), individual_id is required, and the same shared validator
    # is used.
    .validate_individual_timepoint_pairing(
      individual_id = individual_id,
      timepoint     = timepoint,
      n             = n,
      context       = "selection_params$method = 'rise' with rise_paired = TRUE"
    )

    # Sort by individual_id (then timepoint) so that pre-treatment
    # (timepoint == 0) and post-treatment (timepoint == 1) rows for the same
    # individual are matched by position, regardless of input row order.
    rise_order <- order(individual_id, timepoint)
    X                   <- X[rise_order, , drop = FALSE]
    Y                   <- Y[rise_order]
    treatment_pipeline  <- if (!is.null(treatment_pipeline))
      treatment_pipeline[rise_order] else NULL
    covariates_pipeline <- if (!is.null(covariates_pipeline))
      covariates_pipeline[rise_order, , drop = FALSE] else NULL
    individual_id       <- individual_id[rise_order]
    timepoint           <- timepoint[rise_order]

    message(
      "[predictomics] Paired RISE mode detected: rows are sorted by ",
      "individual_id (then timepoint) so that pre-treatment ",
      "(timepoint == 0) and post-treatment (timepoint == 1) rows for the ",
      "same individual are correctly matched, regardless of input row ",
      "order. Predictive modelling will be performed on post-treatment ",
      "samples only (timepoint == 1); pre-treatment samples are used for ",
      "RISE screening only and excluded from the CV loop."
    )
  }

  # ---------------------------------------------------------------------------
  # 1d. selection_params$method = "variance" combined with
  # engineering_params$col_transform = "z" would select on a near-constant
  # (no aggregation) or size/coherence-confounded (with aggregation)
  # post-standardisation variance, since z-scoring fixes each feature's own
  # variance to ~1 before any aggregation. Warn once here (not inside the
  # per-fold loop, which is wrapped in suppressWarnings()) and compute
  # variance scores on the pre-z-score data instead - identical aggregation,
  # if any, just without col_transform - while model fitting continues to
  # use the fully engineered (z-scored) matrix as usual.
  # ---------------------------------------------------------------------------
  is_variance_on_zscore <- !is.null(selection_params) &&
    isTRUE(selection_params$method == "variance") &&
    !is.null(engineering_params) &&
    identical(engineering_params$col_transform %||% "none", "z")

  engineering_params_for_selection <- NULL

  if (is_variance_on_zscore) {
    warning(
      "[predictomics] selection_params$method = 'variance' combined with ",
      "engineering_params$col_transform = 'z' would select features by a ",
      "near-constant or confounded post-standardisation variance (z-scoring ",
      "fixes each feature's own variance to ~1, before any aggregation). ",
      "Variance scores will instead be computed on the pre-z-score data ",
      "(with any geneset aggregation still applied); model fitting continues ",
      "to use the fully engineered (z-scored) matrix as usual.",
      call. = FALSE
    )
    engineering_params_for_selection <- engineering_params_pipeline
    engineering_params_for_selection$col_transform <- "none"
  }

  # ---------------------------------------------------------------------------
  # 2. Outside-CV warning
  # ---------------------------------------------------------------------------
  if (outside_cv) {
    warning(
      "\n[predictomics] outside_cv = TRUE: feature engineering and/or ",
      "selection will be applied to the FULL dataset before cross-validation.\n",
      "  This introduces data leakage and will produce OPTIMISTICALLY BIASED ",
      "performance estimates.\n",
      "  This mode is recommended for exploratory analyses only.\n",
      "  Use outside_cv = FALSE for statistically valid estimates.",
      call. = FALSE, immediate. = TRUE
    )
  }

  # ---------------------------------------------------------------------------
  # 3. Double selection message
  # ---------------------------------------------------------------------------
  embedded_methods <- c("glmnet", "lasso")
  is_embedded      <- model_params$method %in% embedded_methods

  if (!is.null(selection_params) && is_embedded) {
    message(
      "[predictomics] Note: both an explicit selection method ('",
      selection_params$method, "') and an embedded selection model ('",
      model_params$method, "') are specified. Double variable selection will ",
      "be performed: features are first filtered by '", selection_params$method,
      "', then further selected by the regularisation in '",
      model_params$method, "'."
    )
  }

  # ---------------------------------------------------------------------------
  # 3b. compute_importance + scale warning
  # ---------------------------------------------------------------------------
  compute_importance <- model_params$compute_importance %||% FALSE
  coefficient_based_methods <- c("lm", "glmnet", "ridge", "lasso", "svr")

  if (compute_importance && model_params$method %in% coefficient_based_methods &&
      !isTRUE(model_params$scale)) {
    warning(
      "[predictomics] model_params$compute_importance = TRUE with method = '",
      model_params$method, "' but model_params$scale is not TRUE: feature ",
      "importance for this method is a raw (signed) coefficient, which is ",
      "only comparable across features that are on a common scale. Without ",
      "model_params$scale = TRUE, features with larger raw magnitudes will ",
      "tend to receive smaller coefficients (and vice versa), making the ",
      "importance ranking misleading. Set model_params$scale = TRUE for a ",
      "meaningful cross-feature comparison.",
      call. = FALSE, immediate. = TRUE
    )
  }

  # ---------------------------------------------------------------------------
  # 5. Prepare protected predictor matrices (once, outside loop - no leakage)
  # ---------------------------------------------------------------------------
  treatment_mat <- if (!is.null(treatment_pipeline) && treatment_predictor) {
    .prepare_treatment_matrix(treatment_pipeline)
  } else {
    NULL
  }

  covariate_mat <- if (!is.null(covariates_pipeline)) {
    .prepare_covariate_matrix(covariates_pipeline)
  } else {
    NULL
  }

  .validate_predictor_name_collisions(X, treatment_mat, covariate_mat)

  # ---------------------------------------------------------------------------
  # 5b. Generate fold assignments
  # ---------------------------------------------------------------------------
  # Full data (both arms) is retained for passing to run_selection (RISE
  # screening uses both arms). Modelling uses post-treatment only.
  X_processed         <- X   # initialise here for paired RISE subsetting
  X_full              <- X_processed
  Y_full              <- Y
  treatment_full      <- treatment_pipeline
  covariate_mat_full  <- covariate_mat
  treatment_pipeline_post <- treatment_pipeline

  if (is_paired_rise) {
    post_idx      <- which(timepoint == 1)
    X_processed   <- X_processed[post_idx, , drop = FALSE]
    Y             <- Y[post_idx]
    covariate_mat <- if (!is.null(covariate_mat))
      covariate_mat[post_idx, , drop = FALSE]
    else NULL
    treatment_mat <- if (!is.null(treatment_mat))
      treatment_mat[post_idx, , drop = FALSE]
    else NULL
    treatment_pipeline_post <- if (!is.null(treatment_pipeline))
      treatment_pipeline[post_idx] else NULL
    n             <- length(Y)
    folds         <- if (cv_type == "loo") n else folds
  }

  fold_ids <- make_folds(n = n, cv_type = cv_type, k = folds, seed = seed)

  if (verbose) {
    message("[predictomics] Starting ", folds, "-fold CV on ", n,
            " samples and ", p, " features.")
  }

  # ---------------------------------------------------------------------------
  # 6. Outside-CV steps (applied once to full dataset, with leakage)
  # ---------------------------------------------------------------------------
  outside_cv_selection <- NULL

  if (outside_cv) {

    if (!is.null(engineering_params_pipeline)) {
      if (verbose) message("[predictomics] Applying feature engineering outside CV loop.")
      eng_fit  <- run_engineering(
        X_train = if (is_paired_rise) X_full else X,
        params  = engineering_params_pipeline
      )
      if (is_paired_rise) {
        # Store full engineered matrix for RISE selection, then subset for modelling
        X_full      <- eng_fit$X_transformed
        X_processed <- X_full[post_idx, , drop = FALSE]
      } else {
        X_processed <- eng_fit$X_transformed
      }
    }

    if (!is.null(selection_params_pipeline)) {
      if (verbose) message("[predictomics] Applying feature selection outside CV loop.")

      # variance + col_transform = "z": score on the pre-z-score data (X is
      # still raw at this point - only dearseq/gene_level_fc, never
      # col_transform, touch it upstream), same aggregation as the real
      # engineered matrix. X_processed itself is untouched, so modelling
      # still uses the fully engineered (z-scored) matrix below.
      X_for_selection <- if (!is.null(engineering_params_for_selection)) {
        run_engineering(X_train = X, params = engineering_params_for_selection)$X_transformed
      } else if (is_paired_rise) {
        X_full
      } else {
        X_processed
      }

      sel_fit <- run_selection(
        X_train       = X_for_selection,
        Y_train       = if (is_paired_rise) Y_full else Y,
        covariates    = if (is_paired_rise) covariate_mat_full
        else covariate_mat,
        treatment     = if (!is.null(treatment_pipeline))
          .coerce_treatment_binary(
            if (is_paired_rise) treatment_full else treatment_pipeline
          )
        else NULL,
        individual_id = individual_id,
        timepoint     = timepoint,
        params        = selection_params_pipeline
      )
      # Subset selected features to post-treatment X_processed for modelling
      X_processed <- X_processed[, sel_fit$selected_features, drop = FALSE]
      outside_cv_selection <- list(
        selected_features = sel_fit$selected_features,
        selection_scores  = sel_fit$selection_scores,
        n_selected        = length(sel_fit$selected_features)
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 7. Cross-validation loop (parallelised via future_lapply)
  # ---------------------------------------------------------------------------
  if (verbose) message("[predictomics] Running CV folds ...")

  fold_results <- suppressMessages(suppressWarnings(
    future.apply::future_lapply(
      X           = seq_len(folds),
      FUN         = function(k) {
        .run_fold(
          k                  = k,
          fold_ids           = fold_ids,
          X_processed        = X_processed,
          Y                  = Y,
          outside_cv         = outside_cv,
          engineering_params = engineering_params_pipeline,
          engineering_params_for_selection = engineering_params_for_selection,
          selection_params   = selection_params_pipeline,
          model_params       = model_params,
          treatment          = treatment_pipeline,
          treatment_mat      = treatment_mat,
          covariate_mat      = covariate_mat,
          is_paired_rise     = is_paired_rise,
          X_full             = X_full,
          Y_full             = Y_full,
          treatment_full     = treatment_full,
          covariate_mat_full = covariate_mat_full,
          individual_id_full = individual_id,
          timepoint_full     = timepoint
        )
      },
      future.seed = seed
    )
  ))

  # ---------------------------------------------------------------------------
  # 8. Assemble results from fold list
  # ---------------------------------------------------------------------------
  predictions                        <- numeric(n)
  fold_selection_diagnostics         <- if (!is.null(selection_params_pipeline))
    vector("list", folds) else NULL
  fold_embedded_selection_diagnostics <- if (is_embedded)
    vector("list", folds) else NULL
  fold_feature_importance             <- if (compute_importance)
    vector("list", folds) else NULL

  for (k in seq_len(folds)) {
    test_idx              <- which(fold_ids == k)
    predictions[test_idx] <- fold_results[[k]]$predictions[test_idx]

    if (!is.null(selection_params_pipeline))
      fold_selection_diagnostics[[k]] <- fold_results[[k]]$selection_diagnostics

    if (is_embedded)
      fold_embedded_selection_diagnostics[[k]] <-
      fold_results[[k]]$embedded_selection_diagnostics

    if (compute_importance)
      fold_feature_importance[[k]] <- fold_results[[k]]$feature_importance_diagnostics
  }

  if (verbose) message("[predictomics] CV complete.")

  # ---------------------------------------------------------------------------
  # 9. Assemble and return result object
  # ---------------------------------------------------------------------------
  # n_samples/n_samples_total report the pre-restriction (both-arms) count for
  # either paired filtration approach (RISE-paired or dearseq_mode = "paired"),
  # mirroring each other exactly; n_samples_modelled always reports the
  # post-treatment-only count actually used for modelling.
  is_any_paired_filter <- is_paired_rise || is_dearseq_paired
  n_samples_total_value <- if (is_paired_rise) {
    length(Y_full)
  } else if (is_dearseq_paired) {
    n_samples_total_paired
  } else {
    n
  }

  structure(
    list(
      observed                            = Y,
      predicted                           = predictions,
      fold_ids                            = fold_ids,
      treatment                           = if (is_any_paired_filter) treatment_pipeline_post
                                            else treatment,
      treatment_predictor                 = treatment_predictor,
      covariates                          = covariates,
      engineering_params                  = engineering_params,
      selection_params                    = selection_params,
      model_params                        = model_params,
      outside_cv                          = outside_cv,
      cv_type                             = cv_type,
      n_folds                             = folds,
      n_samples                           = n_samples_total_value,
      n_samples_modelled                  = n,
      n_samples_total                     = n_samples_total_value,
      paired_rise                         = is_paired_rise,
      dearseq_paired                      = is_dearseq_paired,
      n_features_input                    = p,
      baseline_model                      = is_baseline_model,
      dearseq_selection                   = dearseq_selection,
      fold_selection_diagnostics          = fold_selection_diagnostics,
      outside_cv_selection                = outside_cv_selection,
      fold_embedded_selection_diagnostics = fold_embedded_selection_diagnostics,
      fold_feature_importance             = fold_feature_importance,
      call                                = cl
    ),
    class = "predictomics"
  )
}
