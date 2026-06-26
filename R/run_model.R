# =============================================================================
# run_model.R
# Model fitting and prediction for the predictomics package.
#
# Implements a unified interface for fitting predictive models with inner-CV
# hyperparameter tuning via the caret package. Supported models:
#   - "lm"     : ordinary least squares (no tuning)
#   - "glmnet" : elastic net regression (alpha, lambda tuned by inner CV)
#   - "ridge"  : ridge regression via glmnet package (lambda tuned by inner CV)
#   - "lasso"  : lasso regression via glmnet package (lambda tuned by inner CV)
#   - "ranger" : random forest via ranger (mtry, min.node.size tuned by inner CV)
#
# run_model()     fits a model on training data, performing inner-CV tuning
#                 where applicable.
# predict_model() applies the fitted model to new data.
# =============================================================================


# -----------------------------------------------------------------------------
#' Fit a predictive model with inner-CV hyperparameter tuning
#'
#' @description
#' Fits a predictive model on the training data \code{X_train} and
#' \code{Y_train}. For models with hyperparameters (\code{"glmnet"},
#' \code{"lasso"}, \code{"ridge"}, \code{"ranger"}), tuning is performed via
#' k-fold cross-validation on the training data using \code{caret::train}, and
#' the best hyperparameters are used to refit a final model. For \code{"lm"},
#' no tuning is performed.
#'
#' @details
#' All models are fitted through \code{caret::train} to provide a uniform
#' interface. Column names of \code{X_train} are sanitised with
#' \code{make.names} prior to fitting to avoid issues with special characters
#' in caret's formula handling. The same sanitisation must be applied to
#' \code{X_new} in \code{\link{predict_model}}, which is handled automatically.
#'
#' Reproducibility of the inner CV is ensured via the \code{future.seed}
#' argument passed to \code{future.apply::future_lapply} in
#' \code{\link{predict_cv}}, which manages per-worker RNG streams. No manual
#' \code{set.seed()} call is made inside \code{run_model} to avoid
#' interfering with these streams.
#'
#' **glmnet/lasso/ridge**: tunes \code{lambda} (regularisation strength) via
#' inner CV. For \code{"glmnet"}, \code{alpha} is also tuned. For
#' \code{"lasso"} and \code{"ridge"}, \code{alpha} is fixed at 1 and 0
#' respectively. If \code{X_train} contains only a single feature, falls back
#' to \code{"lm"} with a warning. For \code{"lasso"} and \code{"glmnet"}
#' (when alpha > 0), non-zero coefficient features are stored in
#' \code{selected_features} of the returned fit object.
#'
#' **ranger**: tunes \code{mtry} (number of variables sampled per split),
#' \code{min.node.size} (minimum node size), and \code{splitrule}
#' (splitting criterion). \code{num.threads} is fixed to 1 to avoid thread
#' oversubscription when the outer CV loop is parallelised.
#'
#' @param X_train Numeric matrix of dimensions n (samples) x p (features).
#'   Training predictor matrix. Must have column names.
#' @param Y_train Numeric vector of length n. Training response variable.
#' @param params A named list of model parameters. See \code{\link{predict_cv}}
#'   for the \code{model_params} argument. Supported fields:
#'   \describe{
#'     \item{\code{method}}{Character string. One of \code{"lm"},
#'       \code{"glmnet"}, \code{"lasso"}, \code{"ridge"}, or \code{"ranger"}.
#'       Required.}
#'     \item{\code{inner_folds}}{Positive integer. Number of inner CV folds
#'       for hyperparameter tuning. Defaults to \code{5}. Ignored for
#'       \code{"lm"}.}
#'     \item{\code{tune_grid}}{A data frame of hyperparameter combinations to
#'       evaluate, or \code{NULL} to use the caret default grid. Optional.}
#'     \item{\code{tune_length}}{Positive integer. Passed to caret's
#'       \code{tuneLength} to auto-generate a grid of this size. Ignored if
#'       \code{tune_grid} is provided. Optional.}
#'     \item{\code{fold_id}}{Integer. Outer fold index, stored in the result
#'       for reference. Defaults to \code{0}. Reproducibility is managed by
#'       \code{future.apply} at the outer CV level.}
#'   }
#'
#' @return A named list of class \code{"predictomics_model"} containing:
#'   \describe{
#'     \item{\code{caret_fit}}{The \code{caret} train object.}
#'     \item{\code{method}}{Character string. The model method as supplied
#'       by the user (e.g. \code{"lasso"}, not the resolved \code{"glmnet"}).}
#'     \item{\code{best_params}}{Data frame. The best hyperparameters selected
#'       by inner CV (\code{NA} for \code{"lm"}).}
#'     \item{\code{inner_folds}}{Integer. Number of inner CV folds used.}
#'     \item{\code{col_names}}{Character vector. Sanitised column names used
#'       during fitting, required for consistent prediction.}
#'     \item{\code{selected_features}}{Character vector of features with
#'       non-zero coefficients at the best lambda, sorted by decreasing
#'       absolute coefficient value. For \code{"lasso"} and \code{"glmnet"}
#'       when alpha > 0. \code{NULL} for all other methods and \code{"ridge"}.}
#'     \item{\code{selection_scores}}{Named numeric vector of raw (signed)
#'       coefficient values for selected features, sorted by decreasing
#'       absolute value. \code{NULL} when \code{selected_features} is
#'       \code{NULL}.}
#'   }
#'
#' @seealso \code{\link{predict_model}}, \code{\link{predict_cv}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X_train <- matrix(rnorm(40 * 100), nrow = 40)
#' colnames(X_train) <- paste0("gene", 1:100)
#' Y_train <- rnorm(40)
#'
#' # OLS
#' fit <- run_model(X_train, Y_train, params = list(method = "lm"))
#'
#' # Lasso with default lambda grid
#' fit <- run_model(X_train, Y_train, params = list(method = "lasso"))
#' fit$selected_features  # non-zero coefficient features
#'
#' # Elastic net with default tuning grid
#' fit <- run_model(X_train, Y_train, params = list(method = "glmnet",
#'                                                   inner_folds = 5))
#'
#' # Random forest with custom tuning grid
#' fit <- run_model(X_train, Y_train, params = list(
#'   method    = "ranger",
#'   tune_grid = expand.grid(mtry = c(5, 10, 20),
#'                           min.node.size = c(1, 5),
#'                           splitrule = "variance")
#' ))
#' }
#'
#' @export
# -----------------------------------------------------------------------------
run_model <- function(X_train, Y_train, params) {

  # ---------------------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------------------
  .validate_X(X_train)
  .validate_Y(Y_train)
  .validate_Y_X_compat(Y_train, X_train)
  .validate_model_params(params, n_train = nrow(X_train))

  if (is.null(colnames(X_train)))
    stop("[predictomics] X_train must have column names for model fitting.",
         call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2. Extract and resolve parameters
  # ---------------------------------------------------------------------------
  method      <- params$method
  inner_folds <- params$inner_folds %||% 5L
  tune_grid   <- params$tune_grid
  tune_length <- params$tune_length %||% 3L
  fold_id     <- params$fold_id     %||% 0L

  # Store the user-supplied method name for the result object
  user_method <- method

  # Resolve lasso/ridge to glmnet with fixed alpha
  fixed_alpha <- NULL
  if (method %in% c("ridge", "lasso")) {
    fixed_alpha <- if (method == "ridge") 0 else 1
    if (is.null(tune_grid))
      tune_grid <- expand.grid(alpha  = fixed_alpha,
                               lambda = 10^seq(-3, 1, length = 50))
    method <- "glmnet"
  }

  # glmnet fallback: requires at least 2 predictors
  if (method == "glmnet" && ncol(X_train) < 2L) {
    warning(
      "[predictomics] glmnet requires at least 2 predictors but X_train has ",
      ncol(X_train), " column(s). Falling back to method = 'lm'.",
      call. = FALSE
    )
    method      <- "lm"
    user_method <- "lm"
  }

  # ---------------------------------------------------------------------------
  # 3. Prepare data: sanitise column names and coerce to data frame
  # ---------------------------------------------------------------------------
  clean_names    <- make.names(colnames(X_train), unique = TRUE)
  X_df           <- as.data.frame(X_train)
  colnames(X_df) <- clean_names
  train_data     <- cbind(X_df, .Y = Y_train)

  # ---------------------------------------------------------------------------
  # 4. Configure caret trainControl
  # ---------------------------------------------------------------------------
  tr_control <- caret::trainControl(
    method          = if (method == "lm") "none" else "cv",
    number          = inner_folds,
    savePredictions = FALSE,
    verboseIter     = FALSE,
    allowParallel   = FALSE
  )

  # ---------------------------------------------------------------------------
  # 5. Fit via caret
  # ---------------------------------------------------------------------------
  caret_fit <- switch(method,

                      lm = {
                        caret::train(
                          .Y ~ .,
                          data      = train_data,
                          method    = "lm",
                          trControl = tr_control
                        )
                      },

                      glmnet = {
                        caret::train(
                          .Y ~ .,
                          data       = train_data,
                          method     = "glmnet",
                          trControl  = tr_control,
                          tuneGrid   = tune_grid,
                          tuneLength = if (is.null(tune_grid)) tune_length else NULL,
                          metric     = "RMSE"
                        )
                      },

                      ranger = {
                        caret::train(
                          .Y ~ .,
                          data        = train_data,
                          method      = "ranger",
                          trControl   = tr_control,
                          tuneGrid    = tune_grid,
                          tuneLength  = if (is.null(tune_grid)) tune_length else NULL,
                          metric      = "RMSE",
                          num.threads = 1L
                        )
                      }
  )

  # ---------------------------------------------------------------------------
  # 6. Extract embedded selected features (lasso / glmnet with alpha > 0)
  # ---------------------------------------------------------------------------
  selected_features <- .extract_glmnet_features(
    caret_fit      = caret_fit,
    method         = method,
    user_method    = user_method,
    fixed_alpha    = fixed_alpha,
    col_names      = clean_names,
    original_names = colnames(X_train)
  )

  # ---------------------------------------------------------------------------
  # 7. Assemble and return fit object
  # ---------------------------------------------------------------------------
  structure(
    list(
      caret_fit         = caret_fit,
      method            = user_method,
      best_params       = if (method == "lm") NA else caret_fit$bestTune,
      inner_folds       = inner_folds,
      col_names         = clean_names,
      selected_features = selected_features$features,
      selection_scores  = selected_features$scores
    ),
    class = "predictomics_model"
  )
}


# -----------------------------------------------------------------------------
#' Generate predictions from a fitted predictomics model
#'
#' @description
#' Applies a model fitted by \code{\link{run_model}} to a new predictor matrix
#' \code{X_new}, returning a numeric vector of predicted values.
#'
#' @param fit A \code{predictomics_model} object returned by
#'   \code{\link{run_model}}.
#' @param X_new Numeric matrix. New predictor matrix. Must have the same
#'   column names as the training matrix passed to \code{\link{run_model}}
#'   (pre-sanitisation).
#'
#' @return A numeric vector of predicted values of length \code{nrow(X_new)}.
#'
#' @seealso \code{\link{run_model}}
#'
#' @export
# -----------------------------------------------------------------------------
predict_model <- function(fit, X_new) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!inherits(fit, "predictomics_model"))
    stop("[predictomics] fit must be a predictomics_model object returned by ",
         "run_model().", call. = FALSE)
  if (!is.matrix(X_new) || !is.numeric(X_new))
    stop("[predictomics] X_new must be a numeric matrix.", call. = FALSE)
  if (is.null(colnames(X_new)))
    stop("[predictomics] X_new must have column names.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2. Sanitise column names to match training
  # ---------------------------------------------------------------------------
  X_df           <- as.data.frame(X_new)
  colnames(X_df) <- make.names(colnames(X_new), unique = TRUE)

  if (!identical(colnames(X_df), fit$col_names))
    stop("[predictomics] Column names of X_new do not match those used during ",
         "model training after sanitisation. Ensure X_new has the same ",
         "features in the same order as X_train.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 3. Predict and return
  # ---------------------------------------------------------------------------
  as.numeric(predict(fit$caret_fit, newdata = X_df))
}


# =============================================================================
# Internal helpers
# =============================================================================

# -----------------------------------------------------------------------------
#' Extract non-zero coefficient features from a fitted glmnet model
#'
#' @description
#' For \code{"lasso"} and \code{"glmnet"} fits where alpha > 0, extracts the
#' names of features with non-zero coefficients at the selected lambda from the
#' final model. Returns \code{NULL} for \code{"ridge"} (alpha = 0), \code{"lm"},
#' and \code{"ranger"}, since these methods do not produce sparse solutions.
#'
#' @param caret_fit A fitted caret train object.
#' @param method Character. The resolved method (\code{"glmnet"} or \code{"lm"}).
#' @param user_method Character. The user-supplied method name.
#' @param fixed_alpha Numeric or \code{NULL}. The fixed alpha value for lasso
#'   or ridge, or \code{NULL} for glmnet (alpha tuned).
#' @param col_names Character vector. Sanitised column names from training.
#' @param original_names Character vector. Original variable names to be passed to output.
#'
#' @return Character vector of selected feature names (in original
#'   pre-sanitisation names where possible), or \code{NULL}.
#' @keywords internal
# -----------------------------------------------------------------------------
.extract_glmnet_features <- function(caret_fit, method, user_method,
                                     fixed_alpha, col_names,
                                     original_names) {

  # Only applicable for glmnet-family with sparsity (alpha > 0)
  is_glmnet  <- method == "glmnet"
  is_ridge   <- !is.null(fixed_alpha) && fixed_alpha == 0
  is_sparse  <- is_glmnet && !is_ridge

  if (!is_sparse) return(list(features = NULL, scores = NULL))

  # Determine best alpha from tuning (for glmnet, alpha may have been tuned)
  best_alpha <- if (!is.null(fixed_alpha)) {
    fixed_alpha
  } else {
    caret_fit$bestTune$alpha
  }

  # Ridge produces no zeros even if alpha = 0 was selected by tuning
  if (!is.null(best_alpha) && best_alpha == 0) return(list(features = NULL, scores = NULL))

  # Extract coefficients at best lambda from the final model
  best_lambda <- caret_fit$bestTune$lambda
  coef_mat    <- coef(caret_fit$finalModel, s = best_lambda)

  # coef_mat is a sparse matrix; convert to named vector and drop intercept
  coef_vec    <- as.numeric(coef_mat)
  coef_names  <- rownames(coef_mat)
  coef_vec    <- setNames(coef_vec, coef_names)
  coef_vec    <- coef_vec[names(coef_vec) != "(Intercept)"]

  # Identify non-zero features using sanitised names (matched to coef_vec)
  # but return original names for interpretability
  nonzero_mask   <- coef_vec != 0
  nonzero_coefs  <- coef_vec[nonzero_mask]
  nonzero_clean  <- col_names[nonzero_mask]
  nonzero_orig   <- original_names[nonzero_mask]

  # Sort by decreasing absolute coefficient value
  ord                 <- order(abs(nonzero_coefs), decreasing = TRUE)
  sorted_orig         <- nonzero_orig[ord]
  sorted_coefs        <- nonzero_coefs[ord]
  names(sorted_coefs) <- sorted_orig

  list(features = sorted_orig, scores = sorted_coefs)
}
