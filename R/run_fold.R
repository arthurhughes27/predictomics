# =============================================================================
# run_fold.R
# Single CV fold execution for the predictomics package.
# Called by future_lapply in predict_cv().
# =============================================================================


# -----------------------------------------------------------------------------
#' Execute a single CV fold
#'
#' @description
#' Internal workhorse called by \code{future.apply::future_lapply} in
#' \code{\link{predict_cv}}. Runs the full pipeline for one fold: engineering,
#' selection, protected predictor appending, model fitting, and prediction.
#' Returns predictions for the test partition, optional explicit selection
#' diagnostics, and optional embedded selection diagnostics.
#'
#' @param k Integer. Fold index.
#' @param fold_ids Integer vector of fold assignments.
#' @param X_processed Numeric matrix. Predictor matrix after any outside-CV
#'   steps.
#' @param Y Numeric vector. Full response vector.
#' @param outside_cv Logical. Whether outside-CV steps have already been
#'   applied.
#' @param engineering_params,selection_params,model_params Pipeline parameter
#'   lists.
#' @param treatment Original treatment vector or \code{NULL}.
#' @param treatment_mat Numeric matrix of treatment columns or \code{NULL}.
#' @param covariate_mat Numeric matrix of covariate columns or \code{NULL}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{predictions}}{Numeric vector of length n with predictions
#'       for test samples in their original positions and zeros elsewhere.}
#'     \item{\code{selection_diagnostics}}{Named list of explicit selection
#'       results for this fold, or \code{NULL}.}
#'     \item{\code{embedded_selection_diagnostics}}{Named list of embedded
#'       selection results (lasso/glmnet non-zero features) for this fold,
#'       or \code{NULL}.}
#'   }
#' @keywords internal
# -----------------------------------------------------------------------------
.run_fold <- function(k, fold_ids, X_processed, Y, outside_cv,
                      engineering_params, selection_params, model_params,
                      treatment, treatment_mat, covariate_mat,
                      is_paired_rise     = FALSE,
                      X_full             = NULL,
                      Y_full             = NULL,
                      treatment_full     = NULL,
                      covariate_mat_full = NULL) {

  n         <- length(Y)
  train_idx <- which(fold_ids != k)
  test_idx  <- which(fold_ids == k)

  X_train <- X_processed[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  X_test  <- X_processed[test_idx,  , drop = FALSE]

  # -- Inside-CV feature engineering ----------------------------------------
  if (!outside_cv && !is.null(engineering_params)) {
    eng_fit <- run_engineering(X_train = X_train, params = engineering_params)
    X_train <- eng_fit$X_transformed
    X_test  <- predict_engineering(eng_fit$fit, X_new = X_test)
  }

  # -- Inside-CV feature selection ------------------------------------------
  selection_diagnostics <- NULL

  if (!outside_cv && !is.null(selection_params)) {
    if (is_paired_rise) {
      # For paired RISE, selection uses the full paired dataset (both arms).
      # The fold indices in the modelling space (post-treatment only) map to
      # post-treatment rows in the full dataset via the stored post_idx.
      # We pass the full training partition of the paired data to run_selection.
      n_full      <- nrow(X_full)
      n_post      <- nrow(X_processed)
      post_idx    <- which(treatment_full == 1)
      pre_idx     <- which(treatment_full == 0)

      # Training post-treatment indices in full data space
      post_train  <- post_idx[train_idx]
      pre_train   <- pre_idx[train_idx]   # matched pre rows for same individuals
      paired_train_idx <- c(pre_train, post_train)

      sel_fit <- run_selection(
        X_train    = X_full[paired_train_idx, , drop = FALSE],
        Y_train    = Y_full[paired_train_idx],
        covariates = if (!is.null(covariate_mat_full))
          covariate_mat_full[paired_train_idx, , drop = FALSE]
        else NULL,
        treatment  = treatment_full[paired_train_idx],
        params     = selection_params
      )
    } else {
      sel_fit <- run_selection(
        X_train    = X_train,
        Y_train    = Y_train,
        covariates = if (!is.null(covariate_mat))
          covariate_mat[train_idx, , drop = FALSE]
        else NULL,
        treatment  = if (!is.null(treatment)) treatment[train_idx] else NULL,
        params     = selection_params
      )
    }
    X_train <- X_train[, sel_fit$selected_features, drop = FALSE]
    X_test  <- X_test[,  sel_fit$selected_features, drop = FALSE]

    selection_diagnostics <- list(
      selected_features = sel_fit$selected_features,
      selection_scores  = sel_fit$selection_scores,
      n_selected        = length(sel_fit$selected_features)
    )
  }

  # -- Append protected predictors (treatment then covariates) ---------------
  if (!is.null(treatment_mat)) {
    X_train <- cbind(X_train, treatment_mat[train_idx, , drop = FALSE])
    X_test  <- cbind(X_test,  treatment_mat[test_idx,  , drop = FALSE])
  }

  if (!is.null(covariate_mat)) {
    X_train <- cbind(X_train, covariate_mat[train_idx, , drop = FALSE])
    X_test  <- cbind(X_test,  covariate_mat[test_idx,  , drop = FALSE])
  }

  # -- Model fitting and prediction -----------------------------------------
  model_fit <- run_model(
    X_train = X_train,
    Y_train = Y_train,
    params  = modifyList(model_params, list(fold_id = k))
  )

  # -- Embedded selection diagnostics ---------------------------------------
  embedded_selection_diagnostics <- if (!is.null(model_fit$selected_features)) {
    list(
      selected_features = model_fit$selected_features,
      selection_scores  = model_fit$selection_scores,
      n_selected        = length(model_fit$selected_features)
    )
  } else {
    NULL
  }

  # -- Assemble predictions in original sample order ------------------------
  predictions           <- numeric(n)
  predictions[test_idx] <- predict_model(fit = model_fit, X_new = X_test)

  list(
    predictions                    = predictions,
    selection_diagnostics          = selection_diagnostics,
    embedded_selection_diagnostics = embedded_selection_diagnostics
  )
}
