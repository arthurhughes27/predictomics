# =============================================================================
# make_folds.R
# Fold assignment generation for cross-validation in the predictomics package.
# =============================================================================


# -----------------------------------------------------------------------------
#' Generate cross-validation fold assignments
#'
#' @description
#' Generates a vector of fold assignments for use in the cross-validation loop
#' of \code{\link{predict_cv}}. Supports k-fold and leave-one-out (LOO) CV.
#'
#' @details
#' For k-fold CV, samples are randomly assigned to one of \code{k} folds of
#' approximately equal size. The assignment is reproducible given the same
#' \code{seed}.
#'
#' For LOO CV, each sample is assigned its own fold (fold ID equal to its
#' index), so that each iteration of the CV loop holds out exactly one sample.
#' The \code{seed} and \code{k} arguments are ignored in this case.
#'
#' The returned vector of fold IDs is used in \code{\link{predict_cv}} to
#' partition samples into training and test sets: for fold \code{k}, samples
#' with \code{fold_ids == k} form the test set and the remainder form the
#' training set.
#'
#' @param n Positive integer. Number of samples.
#' @param cv_type Character string. Type of cross-validation. One of
#'   \code{"kfold"} (K-fold CV) or \code{"loo"} (leave-one-out CV). Defaults to \code{"kfold"}.
#' @param k Positive integer. Number of folds for k-fold CV. Must satisfy
#'   \code{2 <= k <= n}. Ignored when \code{cv_type = "loo"}.
#'   Defaults to \code{10}.
#' @param seed Integer. Random seed for reproducible fold assignment in k-fold
#'   CV. Ignored when \code{cv_type = "loo"}. Defaults to \code{12345}.
#'
#' @return An integer vector of length \code{n} containing fold assignments.
#'   Values range from \code{1} to \code{k} for k-fold CV, or \code{1} to
#'   \code{n} for LOO CV.
#'
#' @seealso \code{\link{predict_cv}}
#'
#' @examples
#' # K-fold CV: 50 samples into 10 folds
#' fold_ids <- make_folds(n = 50, cv_type = "kfold", k = 10, seed = 12345)
#' table(fold_ids)  # folds should be approximately equal in size
#'
#' # LOO CV: each sample is its own fold
#' fold_ids <- make_folds(n = 20, cv_type = "loo")
#' stopifnot(identical(fold_ids, 1:20))
#'
#' @export
# -----------------------------------------------------------------------------
make_folds <- function(n,
                       cv_type = "kfold",
                       k       = 10L,
                       seed    = 12345L) {

  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  .validate_scalar_args(cv_type = cv_type, folds = k, n = n,
                        seed = seed, outside_cv = FALSE, verbose = TRUE)

  # ---------------------------------------------------------------------------
  # 2. Generate fold assignments
  # ---------------------------------------------------------------------------
  fold_ids <- switch(cv_type,

    kfold = {
      set.seed(seed)
      sample(rep(seq_len(k), length.out = n))
    },

    loo = seq_len(n)
  )

  as.integer(fold_ids)
}
