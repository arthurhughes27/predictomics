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

  cov_df <- as.data.frame(covariates)

  # Detect and remove null covariates (single-level variables)
  n_levels <- vapply(cov_df, function(col) {
    if (is.numeric(col)) length(unique(col[!is.na(col)]))
    else                  nlevels(factor(col[!is.na(col)]))
  }, integer(1))

  null_covs <- names(n_levels)[n_levels < 2L]

  if (length(null_covs) > 0L) {
    warning(
      "[predictomics] The following covariate(s) have only one unique value ",
      "and will be removed from modelling: ",
      paste(null_covs, collapse = ", "), ".",
      call. = FALSE
    )
    cov_df <- cov_df[, n_levels >= 2L, drop = FALSE]
  }

  if (ncol(cov_df) == 0L) {
    message(
      "[predictomics] All covariates were removed due to having only one ",
      "unique value. Proceeding without covariates."
    )
    return(NULL)
  }

  mat <- model.matrix(~ ., data = cov_df)[, -1L, drop = FALSE]

  if (!is.null(rownames(covariates))) rownames(mat) <- rownames(covariates)

  mat
}


# -----------------------------------------------------------------------------
#' Coerce a treatment variable to a binary numeric vector (0/1)
#'
#' @description
#' Safely converts a treatment variable to a binary numeric vector regardless
#' of whether it was supplied as a binary numeric or a two-level factor.
#' For binary numeric input, values are returned unchanged. For a factor,
#' the first level is mapped to 0 and the second level to 1, matching the
#' convention used by \code{model.matrix} and \code{.prepare_treatment_matrix}.
#'
#' This function is used in the outside-CV branch of \code{\link{predict_cv}}
#' before passing treatment to \code{\link{run_selection}} for the
#' \code{"rise"} method, which requires a strict 0/1 numeric vector.
#'
#' @param treatment A binary numeric vector (0/1) or a two-level factor.
#' @return A binary integer vector with values 0 and 1.
#' @keywords internal
# -----------------------------------------------------------------------------
.coerce_treatment_binary <- function(treatment) {
  if (is.factor(treatment)) {
    as.integer(treatment) - 1L
  } else {
    as.integer(treatment)
  }
}
