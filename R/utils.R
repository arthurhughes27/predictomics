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
