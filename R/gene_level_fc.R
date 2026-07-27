# =============================================================================
# gene_level_fc.R
# Gene-level fold-change (post - pre) engineering for paired designs.
# =============================================================================


# -----------------------------------------------------------------------------
#' Collapse paired pre/post observations into individual-level fold changes
#'
#' @description
#' Internal workhorse for \code{gene_level_fc} engineering, called from
#' \code{\link{predict_cv}} before the CV loop and before any other
#' engineering or selection step. Assumes
#' \code{\link{.validate_gene_level_fc_inputs}} has already been run, so every
#' individual has exactly one pre-treatment (\code{timepoint = 0}) and one
#' post-treatment (\code{timepoint = 1}) row.
#'
#' @param X Numeric matrix of dimensions n x p.
#' @param Y Numeric vector of length n.
#' @param individual_id Vector of length n identifying individuals.
#' @param timepoint Binary numeric vector of length n (0/1).
#' @param treatment Optional vector/factor of length n (e.g. the
#'   \code{treatment} argument of \code{\link{predict_cv}}). If supplied, it
#'   is reduced to one entry per individual by taking the post-treatment
#'   (\code{timepoint = 1}) row. \code{NULL} (default) if not supplied.
#' @param covariates Optional matrix or data frame with n rows (e.g. the
#'   \code{covariates} argument of \code{\link{predict_cv}}). If supplied, it
#'   is reduced to one row per individual by taking the post-treatment
#'   (\code{timepoint = 1}) row. \code{NULL} (default) if not supplied.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{X}}{Numeric matrix of dimensions n_individuals x p, each
#'       row equal to the post-treatment minus pre-treatment value for that
#'       individual.}
#'     \item{\code{Y}}{Numeric vector of length n_individuals, taken from the
#'       post-treatment observation of each individual.}
#'     \item{\code{treatment}}{\code{treatment} reduced to one entry per
#'       individual (post-treatment row), or \code{NULL} if not supplied.}
#'     \item{\code{covariates}}{\code{covariates} reduced to one row per
#'       individual (post-treatment row), or \code{NULL} if not supplied.}
#'   }
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_gene_level_fc <- function(X, Y, individual_id, timepoint,
                                   treatment = NULL, covariates = NULL) {

  uid   <- unique(individual_id)
  n_ind <- length(uid)

  X_fc      <- matrix(NA_real_, nrow = n_ind, ncol = ncol(X),
                      dimnames = list(NULL, colnames(X)))
  Y_fc      <- numeric(n_ind)
  post_rows <- integer(n_ind)

  for (i in seq_len(n_ind)) {
    idx      <- which(individual_id == uid[i])
    pre_idx  <- idx[timepoint[idx] == 0]
    post_idx <- idx[timepoint[idx] == 1]

    X_fc[i, ]    <- X[post_idx, ] - X[pre_idx, ]
    Y_fc[i]      <- Y[post_idx]
    post_rows[i] <- post_idx
  }

  rownames(X_fc) <- as.character(uid)

  list(
    X          = X_fc,
    Y          = Y_fc,
    treatment  = if (!is.null(treatment)) treatment[post_rows] else NULL,
    covariates = if (!is.null(covariates)) covariates[post_rows, , drop = FALSE] else NULL
  )
}
