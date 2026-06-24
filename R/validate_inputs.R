# =============================================================================
# validate_inputs.R
# Centralised input validation helpers for the predictomics package.
#
# Design principle: separate validators are defined for arguments with
# non-trivial logic or that are reused across multiple functions. Simple
# scalar flag checks are grouped into .validate_scalar_args() to avoid
# excessive fragmentation.
#
# To extend validation when new arguments are added, either add a new
# validator for complex arguments or extend .validate_scalar_args() for
# simple flags. Do not duplicate validation logic in individual function
# scripts.
# =============================================================================


# -----------------------------------------------------------------------------
#' Validate the response vector Y
#'
#' @description
#' Checks that \code{Y} is a non-empty numeric vector with no \code{NA} values.
#'
#' @param Y The response variable passed to a predictomics function.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_Y <- function(Y) {

  if (!is.numeric(Y) || !is.null(dim(Y)))
    stop("[predictomics] Y must be a numeric vector.", call. = FALSE)
  if (length(Y) < 2L)
    stop("[predictomics] Y must contain at least 2 observations.", call. = FALSE)
  if (anyNA(Y))
    stop("[predictomics] Y contains NA values. Please impute or remove them.",
         call. = FALSE)

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Validate the predictor matrix X
#'
#' @description
#' Checks that \code{X} is a numeric matrix with at least one feature and no
#' \code{NA} values.
#'
#' @param X The predictor matrix passed to a predictomics function.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_X <- function(X) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("[predictomics] X must be a numeric matrix.", call. = FALSE)
  if (nrow(X) < 2L || ncol(X) < 1L)
    stop("[predictomics] X must have at least 2 rows (samples) and 1 column ",
         "(feature).", call. = FALSE)
  if (anyNA(X))
    stop("[predictomics] X contains NA values. Please impute or remove them.",
         call. = FALSE)

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Validate compatibility between Y and X
#'
#' @description
#' Checks that the number of observations in \code{Y} matches the number of
#' rows in \code{X}. Call after \code{.validate_Y} and \code{.validate_X}.
#'
#' @param Y The response variable passed to a predictomics function.
#' @param X The predictor matrix passed to a predictomics function.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_Y_X_compat <- function(Y, X) {

  if (length(Y) != nrow(X))
    stop("[predictomics] length(Y) (", length(Y), ") must equal nrow(X) (",
         nrow(X), ").", call. = FALSE)

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Validate a pipeline parameter list
#'
#' @description
#' Checks that a parameter list is either \code{NULL} or a named list
#' containing a \code{method} character element. Used for
#' \code{engineering_params}, \code{selection_params}, and
#' \code{model_params}.
#'
#' @param params The parameter list to validate.
#' @param arg_name Character string. The argument name, used in error messages.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_params_list <- function(params, arg_name) {

  if (is.null(params)) return(invisible(NULL))

  if (!is.list(params))
    stop("[predictomics] ", arg_name, " must be a named list or NULL.",
         call. = FALSE)
  if (!is.character(params$method) || length(params$method) != 1L)
    stop("[predictomics] ", arg_name, " must contain a 'method' element as a ",
         "single character string.", call. = FALSE)

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Validate simple scalar arguments
#'
#' @description
#' Checks \code{folds}, \code{seed}, \code{outside_cv}, and \code{verbose}.
#' These are grouped together as their checks are straightforward one-liners
#' that do not warrant individual functions.
#'
#' @param cv_type type of CV ("kfold" or "loo")
#' @param folds Number of CV folds.
#' @param n Integer. Number of samples, used to bound \code{folds}.
#' @param seed Random seed.
#' @param outside_cv Logical flag.
#' @param verbose Logical flag.
#' @param treatment_predictor Logical flag.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_scalar_args <- function(cv_type, folds, n, seed, outside_cv, verbose,
                                  treatment_predictor = FALSE) {

  if (!cv_type %in% c("kfold", "loo")){
    stop("[predictomics] cv_type must be one of `kfold` (K-fold CV) or `loo` (leave-one-out CV)",
         call. = FALSE)
  }

  if (!is.numeric(folds) || length(folds) != 1L || folds != as.integer(folds)
      || folds < 2L || folds > n)
    stop("[predictomics] folds must be an integer >= 2 and <= n (", n, ").",
         call. = FALSE)

  if (!is.numeric(seed) || length(seed) != 1L)
    stop("[predictomics] seed must be a single numeric value.", call. = FALSE)

  if (!is.logical(outside_cv) || length(outside_cv) != 1L)
    stop("[predictomics] outside_cv must be TRUE or FALSE.", call. = FALSE)

  if (!is.logical(verbose) || length(verbose) != 1L)
    stop("[predictomics] verbose must be TRUE or FALSE.", call. = FALSE)

  if (!is.logical(treatment_predictor) || length(treatment_predictor) != 1L)
    stop("[predictomics] treatment_predictor must be TRUE or FALSE.",
         call. = FALSE)

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Validate the treatment variable
#'
#' @description
#' Checks that \code{treatment} is either \code{NULL}, a factor with at least
#' 2 levels, or a binary numeric vector containing only 0 and 1. Also checks
#' length compatibility with \code{Y} and absence of \code{NA} values.
#'
#' @param treatment The treatment variable passed to a predictomics function.
#' @param Y The response vector, used for length compatibility check.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_treatment <- function(treatment, Y) {

  if (is.null(treatment)) return(invisible(NULL))

  if (length(treatment) != length(Y))
    stop("[predictomics] treatment must have the same length as Y (",
         length(Y), ").", call. = FALSE)

  if (anyNA(treatment))
    stop("[predictomics] treatment contains NA values.", call. = FALSE)

  if (is.factor(treatment)) {
    if (nlevels(treatment) < 2L)
      stop("[predictomics] treatment must have at least 2 levels.", call. = FALSE)
  } else if (is.numeric(treatment)) {
    if (!all(treatment %in% c(0, 1)))
      stop("[predictomics] Numeric treatment must contain only 0 and 1.",
           call. = FALSE)
  } else {
    stop("[predictomics] treatment must be a factor or a binary numeric ",
         "vector (0/1).", call. = FALSE)
  }

  invisible(NULL)
}



# -----------------------------------------------------------------------------
#' Validate model_params
#'
#' @description
#' Checks that \code{model_params} specifies a supported method and that
#' tuning-related arguments are consistent and valid.
#'
#' @param params The \code{model_params} list.
#' @param n_train Integer. Number of training samples, used to bound
#'   \code{inner_folds}.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_model_params <- function(params, n_train) {

  supported <- c("lm", "glmnet", "ranger")
  if (!params$method %in% supported)
    stop("[predictomics] model_params$method must be one of: ",
         paste(supported, collapse = ", "), ".", call. = FALSE)

  inner_folds <- params$inner_folds %||% 5L
  if (!is.numeric(inner_folds) || length(inner_folds) != 1L ||
      inner_folds != as.integer(inner_folds) || inner_folds < 2L)
    stop("[predictomics] model_params$inner_folds must be an integer >= 2.",
         call. = FALSE)
  if (inner_folds >= n_train)
    stop("[predictomics] model_params$inner_folds (", inner_folds, ") must be ",
         "less than the number of training samples (", n_train, ").",
         call. = FALSE)

  if (!is.null(params$tune_grid) && !is.data.frame(params$tune_grid))
    stop("[predictomics] model_params$tune_grid must be a data frame or NULL.",
         call. = FALSE)

  tune_length <- params$tune_length
  if (!is.null(tune_length)) {
    if (!is.numeric(tune_length) || length(tune_length) != 1L ||
        tune_length != as.integer(tune_length) || tune_length < 1L)
      stop("[predictomics] model_params$tune_length must be a positive integer.",
           call. = FALSE)
  }

  seed    <- params$seed    %||% 12345L
  fold_id <- params$fold_id %||% 0L
  if (!is.numeric(seed)    || length(seed)    != 1L)
    stop("[predictomics] model_params$seed must be a single numeric value.",
         call. = FALSE)
  if (!is.numeric(fold_id) || length(fold_id) != 1L)
    stop("[predictomics] model_params$fold_id must be a single numeric value.",
         call. = FALSE)

  invisible(NULL)
}


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


# -----------------------------------------------------------------------------
#' Validate the covariates matrix or data frame
#'
#' @description
#' Checks that \code{covariates} is either \code{NULL}, a numeric matrix, or
#' a data frame. Verifies row count matches \code{Y}, that column names are
#' present, and that no column is entirely \code{NA}.
#'
#' @param covariates The covariates argument passed to a predictomics function.
#' @param Y The response vector, used for length compatibility check.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_covariates <- function(covariates, Y) {

  if (is.null(covariates)) return(invisible(NULL))

  if (!is.matrix(covariates) && !is.data.frame(covariates))
    stop("[predictomics] covariates must be a numeric matrix or data frame.",
         call. = FALSE)

  if (nrow(covariates) != length(Y))
    stop("[predictomics] covariates must have the same number of rows as ",
         "length(Y) (", length(Y), ").", call. = FALSE)

  if (is.null(colnames(covariates)))
    stop("[predictomics] covariates must have column names.", call. = FALSE)

  all_na <- vapply(seq_len(ncol(covariates)),
                   function(j) all(is.na(covariates[, j])),
                   logical(1))
  if (any(all_na))
    stop("[predictomics] The following covariate column(s) are entirely NA: ",
         paste(colnames(covariates)[all_na], collapse = ", "), ".",
         call. = FALSE)

  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Validate for column name collisions across predictor matrices
#'
#' @description
#' Checks that column names of \code{X}, \code{treatment_mat}, and
#' \code{covariate_mat} do not overlap, to prevent ambiguity when they are
#' column-bound before model fitting.
#'
#' @param X The (possibly engineered/selected) predictor matrix.
#' @param treatment_mat Numeric matrix of treatment columns, or \code{NULL}.
#' @param covariate_mat Numeric matrix of covariate columns, or \code{NULL}.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_predictor_name_collisions <- function(X, treatment_mat,
                                                covariate_mat) {

  x_names   <- colnames(X)
  trt_names <- if (!is.null(treatment_mat)) colnames(treatment_mat) else character(0)
  cov_names <- if (!is.null(covariate_mat)) colnames(covariate_mat) else character(0)

  x_trt <- intersect(x_names, trt_names)
  if (length(x_trt) > 0L)
    stop("[predictomics] Column name collision between X and treatment: ",
         paste(x_trt, collapse = ", "), ".", call. = FALSE)

  x_cov <- intersect(x_names, cov_names)
  if (length(x_cov) > 0L)
    stop("[predictomics] Column name collision between X and covariates: ",
         paste(x_cov, collapse = ", "), ".", call. = FALSE)

  trt_cov <- intersect(trt_names, cov_names)
  if (length(trt_cov) > 0L)
    stop("[predictomics] Column name collision between treatment and ",
         "covariates: ", paste(trt_cov, collapse = ", "), ".", call. = FALSE)

  invisible(NULL)
}
