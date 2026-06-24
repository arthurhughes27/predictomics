# =============================================================================
# utils.R
# Shared utility functions for the predictomics package.
#
# Functions defined here are available to all other scripts in the package.
# Do not export these unless there is a specific reason for users to call them
# directly.
# =============================================================================


# -----------------------------------------------------------------------------
#' Null coalescing operator
#'
#' @name null_coalesce
#' @rdname null_coalesce
#'
#' @description
#' Returns \code{rhs} if \code{lhs} is \code{NULL}, otherwise returns
#' \code{lhs}. Used throughout the package to resolve default parameter values
#' from user-supplied lists cleanly.
#'
#' @param lhs The value to test.
#' @param rhs The fallback value if \code{lhs} is \code{NULL}.
#'
#' @return \code{lhs} if not \code{NULL}, otherwise \code{rhs}.
#'
#' @examples
#' \dontrun{
#' NULL %||% 5     # returns 5
#' 3    %||% 5     # returns 3
#' }
#'
#' @keywords internal
# -----------------------------------------------------------------------------
`%||%` <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs


# -----------------------------------------------------------------------------
#' Prepare a numeric treatment matrix for use as model predictors
#'
#' @description
#' Converts a treatment variable into a numeric matrix suitable for
#' column-binding onto a predictor matrix immediately before model fitting.
#'
#' @details
#' Two cases are handled:
#' \itemize{
#'   \item \strong{Binary numeric} (\code{0}/\code{1}): returned as a single
#'     column matrix named \code{"treatment"}. No encoding is required.
#'   \item \strong{Factor with k levels}: encoded into k-1 dummy columns via
#'     \code{model.matrix(~ treatment)}, which drops the first level as the
#'     reference category. Column names take the form
#'     \code{"treatment<level>"}, e.g. \code{"treatmentActive"}.
#' }
#' This function is called once before the CV loop in \code{\link{predict_cv}}
#' and the resulting matrix is subset by fold index inside the loop. Because
#' treatment assignment is known for all samples (including test fold samples)
#' and is not estimated from data, there is no data leakage from computing this
#' matrix on the full dataset.
#'
#' @param treatment A factor or binary numeric vector of length n.
#'
#' @return A numeric matrix of dimensions n x (1 or k-1) with informative
#'   column names. Row names match those of \code{treatment} if present.
#'
#' @seealso \code{\link{predict_cv}}
#'
#' @keywords internal
# -----------------------------------------------------------------------------
.prepare_treatment_matrix <- function(treatment) {

  if (is.numeric(treatment)) {

    # Binary numeric: single column, no encoding needed
    mat        <- matrix(treatment, ncol = 1L,
                         dimnames = list(names(treatment), "treatment"))

  } else {

    # Factor: k-1 dummy columns, reference level dropped automatically
    # model.matrix returns an intercept column which we discard
    mat <- model.matrix(~ treatment)[, -1L, drop = FALSE]

    # Clean up row names (model.matrix uses seq_len(n) by default)
    if (!is.null(names(treatment))) rownames(mat) <- names(treatment)
  }

  mat
}


# -----------------------------------------------------------------------------
#' Prepare a numeric covariate matrix for use as model predictors
#'
#' @description
#' Converts a covariate matrix or data frame into a numeric matrix suitable
#' for column-binding onto the predictor matrix immediately before model
#' fitting. Factor and character columns are one-hot encoded via
#' \code{model.matrix}, with the first level of each factor dropped as the
#' reference category. Numeric columns are passed through unchanged.
#'
#' @details
#' Encoding is performed once on the full covariate matrix before the CV loop.
#' The resulting numeric matrix is then subset by fold index inside the loop.
#' Because covariate values are fixed design variables (not estimated from
#' data), computing the encoding on the full dataset introduces no leakage.
#'
#' Column names of the returned matrix take the form produced by
#' \code{model.matrix}: numeric columns retain their original name; factor
#' columns produce names of the form \code{"<column><level>"}.
#'
#' @param covariates A numeric matrix or data frame of dimensions n x q.
#'   Must have column names. Factor and character columns are one-hot encoded.
#'
#' @return A numeric matrix of dimensions n x q' (where q' >= q due to
#'   one-hot encoding of factors), with informative column names.
#'
#' @seealso \code{\link{predict_cv}}
#'
#' @keywords internal
# -----------------------------------------------------------------------------
.prepare_covariate_matrix <- function(covariates) {

  # model.matrix handles both numeric and factor columns correctly.
  # The intercept column (first column) is dropped via [, -1L].
  # ~ . expands all columns in the data frame.
  mat <- model.matrix(~ ., data = as.data.frame(covariates))[, -1L, drop = FALSE]

  # Restore original row names if present
  if (!is.null(rownames(covariates))) rownames(mat) <- rownames(covariates)

  mat
}
