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
  if (!is.numeric(Y) || !is.null(dim(Y))) {
    stop("[predictomics] Y must be a numeric vector.", call. = FALSE)
  }
  if (length(Y) < 2L) {
    stop("[predictomics] Y must contain at least 2 observations.",
         call. = FALSE)
  }
  if (anyNA(Y)) {
    stop("[predictomics] Y contains NA values. Please impute or remove them.",
         call. = FALSE)
  }

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
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("[predictomics] X must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 1L) {
    stop(
      "[predictomics] X must have at least 2 rows (samples) and 1 column ",
      "(feature).",
      call. = FALSE
    )
  }
  if (anyNA(X)) {
    stop("[predictomics] X contains NA values. Please impute or remove them.",
         call. = FALSE)
  }

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
  if (length(Y) != nrow(X)) {
    stop(
      "[predictomics] length(Y) (",
      length(Y),
      ") must equal nrow(X) (",
      nrow(X),
      ").",
      call. = FALSE
    )
  }

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
  if (is.null(params))
    return(invisible(NULL))

  if (!is.list(params)) {
    stop("[predictomics] ",
         arg_name,
         " must be a named list or NULL.",
         call. = FALSE)
  }

  if (!is.character(params$method) || length(params$method) != 1L) {
    stop(
      "[predictomics] ",
      arg_name,
      " must contain a 'method' element as a ",
      "single character string.",
      call. = FALSE
    )
  }

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
#' @param folds Number of CV folds.
#' @param n Integer. Number of samples, used to bound \code{folds}.
#' @param seed Random seed.
#' @param outside_cv Logical flag.
#' @param verbose Logical flag.
#' @return Invisibly returns \code{NULL} if validation passes.
#' @keywords internal
# -----------------------------------------------------------------------------
.validate_scalar_args <- function(folds,
                                  n,
                                  seed,
                                  cv_type,
                                  outside_cv,
                                  verbose) {
  if (!is.numeric(folds) ||
      length(folds) != 1L || folds != as.integer(folds)
      || folds < 2L || folds > n) {
    stop("[predictomics] folds must be an integer >= 2 and <= n (",
         n,
         ").",
         call. = FALSE)
  }

  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("[predictomics] seed must be a single numeric value.", call. = FALSE)
  }

  if (!is.logical(outside_cv) || length(outside_cv) != 1L) {
    stop("[predictomics] outside_cv must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("[predictomics] verbose must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.character(cv_type) || !(cv_type %in% c("kfold", "loo"))) {
    stop("[predictomics] cv_type must be one of `kfold` or `loo`.",
         call. = FALSE)
  }

  invisible(NULL)
}
