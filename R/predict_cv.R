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
#'   The predictor matrix. Column names should be feature identifiers.
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
#'   options. See \code{\link{run_model}} for supported options. Defaults to
#'   \code{list(method = "lm")} (ordinary least squares).
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
#'       \code{NULL}.}
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
#'     \item{\code{n_samples}}{Integer. Number of samples.}
#'     \item{\code{n_features_input}}{Integer. Number of features in the input
#'       \code{X}.}
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
#' }
#'
#' @export
# -----------------------------------------------------------------------------
predict_cv <- function(Y,
                       X,
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
                       verbose             = TRUE) {

  cl <- match.call()

  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  n <- length(Y)
  p <- ncol(X)
  folds <- if (cv_type == "loo") n else folds

  .validate_Y(Y)
  .validate_X(X)
  .validate_Y_X_compat(Y, X)
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
  # 4. Generate fold assignments
  # ---------------------------------------------------------------------------
  fold_ids <- make_folds(n = n, cv_type = cv_type, k = folds, seed = seed)

  if (verbose) {
    message("[predictomics] Starting ", folds, "-fold CV on ", n,
            " samples and ", p, " features.")
  }

  # ---------------------------------------------------------------------------
  # 5. Prepare protected predictor matrices (once, outside loop — no leakage)
  # ---------------------------------------------------------------------------
  treatment_mat <- if (!is.null(treatment) && treatment_predictor) {
    .prepare_treatment_matrix(treatment)
  } else {
    NULL
  }

  covariate_mat <- if (!is.null(covariates)) {
    .prepare_covariate_matrix(covariates)
  } else {
    NULL
  }

  .validate_predictor_name_collisions(X, treatment_mat, covariate_mat)

  # ---------------------------------------------------------------------------
  # 6. Outside-CV steps (applied once to full dataset, with leakage)
  # ---------------------------------------------------------------------------
  outside_cv_selection <- NULL
  X_processed          <- X

  if (outside_cv) {

    if (!is.null(engineering_params)) {
      if (verbose) message("[predictomics] Applying feature engineering outside CV loop.")
      eng_fit     <- run_engineering(X_train = X, params = engineering_params)
      X_processed <- eng_fit$X_transformed
    }

    if (!is.null(selection_params)) {
      if (verbose) message("[predictomics] Applying feature selection outside CV loop.")
      sel_fit     <- run_selection(
        X_train    = X_processed,
        Y_train    = Y,
        covariates = covariate_mat,
        treatment  = if (!is.null(treatment))
          .coerce_treatment_binary(treatment)
        else NULL,
        params     = selection_params
      )
      X_processed          <- X_processed[, sel_fit$selected_features, drop = FALSE]
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
          engineering_params = engineering_params,
          selection_params   = selection_params,
          model_params       = model_params,
          treatment          = treatment,
          treatment_mat      = treatment_mat,
          covariate_mat      = covariate_mat
        )
      },
      future.seed = seed
    )
  ))

  # ---------------------------------------------------------------------------
  # 8. Assemble results from fold list
  # ---------------------------------------------------------------------------
  predictions                        <- numeric(n)
  fold_selection_diagnostics         <- if (!is.null(selection_params))
    vector("list", folds) else NULL
  fold_embedded_selection_diagnostics <- if (is_embedded)
    vector("list", folds) else NULL

  for (k in seq_len(folds)) {
    test_idx              <- which(fold_ids == k)
    predictions[test_idx] <- fold_results[[k]]$predictions[test_idx]

    if (!is.null(selection_params))
      fold_selection_diagnostics[[k]] <- fold_results[[k]]$selection_diagnostics

    if (is_embedded)
      fold_embedded_selection_diagnostics[[k]] <-
      fold_results[[k]]$embedded_selection_diagnostics
  }

  if (verbose) message("[predictomics] CV complete.")

  # ---------------------------------------------------------------------------
  # 9. Assemble and return result object
  # ---------------------------------------------------------------------------
  structure(
    list(
      observed                            = Y,
      predicted                           = predictions,
      fold_ids                            = fold_ids,
      treatment                           = treatment,
      treatment_predictor                 = treatment_predictor,
      covariates                          = covariates,
      engineering_params                  = engineering_params,
      selection_params                    = selection_params,
      model_params                        = model_params,
      outside_cv                          = outside_cv,
      cv_type                             = cv_type,
      n_folds                             = folds,
      n_samples                           = n,
      n_features_input                    = p,
      fold_selection_diagnostics          = fold_selection_diagnostics,
      outside_cv_selection                = outside_cv_selection,
      fold_embedded_selection_diagnostics = fold_embedded_selection_diagnostics,
      call                                = cl
    ),
    class = "predictomics"
  )
}
