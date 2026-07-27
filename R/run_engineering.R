# =============================================================================
# run_engineering.R
# Feature engineering for the predictomics package.
#
# Implements three sequential steps:
#   0. Fold-change (FC) engineering: computes per-individual differences
#      between post-treatment and pre-treatment observations. Requires
#      treatment and participant_id vectors. Reduces n_rows by half.
#   1. Column-wise transformation of individual features (z-score or none).
#   2. Geneset aggregation: summarising groups of features into single values
#      (mean, median, sum, or pc1).
#
# Step 0 always precedes steps 1 and 2. Steps 1 and 2 store the parameters
# fitted on the training set so that they can be applied to the test set
# without data leakage via predict_engineering().
# =============================================================================


# -----------------------------------------------------------------------------
#' Apply feature engineering to a training matrix
#'
#' @description
#' Applies a sequential feature engineering pipeline to the predictor matrix
#' \code{X_train}: (0) optional fold-change engineering collapsing pre/post
#' pairs into difference scores, (1) optional column-wise transformation of
#' individual features, followed by (2) optional aggregation of features into
#' genesets. All parameters required to apply the same transformations to a
#' test set are stored in the returned fit object and consumed by
#' \code{\link{predict_engineering}}.
#'
#' @details
#' **Step 0 - Fold-change engineering** collapses paired pre/post observations
#' into a single difference score per individual per feature:
#' \eqn{FC_{ij} = X_{post,ij} - X_{pre,ij}}. Requires \code{treatment} and
#' \code{participant_id} to be supplied. The output has \code{n_individuals}
#' rows (half of the input rows). This step is parameter-free: no training
#' statistics are stored, and the same operation is applied identically to
#' test data.
#'
#' **Step 1 - Column-wise transformation** is applied independently to each
#' feature (column):
#' \itemize{
#'   \item \code{"none"}: no transformation is applied.
#'   \item \code{"z"}: each feature is z-scored using the mean and standard
#'     deviation computed from \code{X_train}. Features with zero variance in
#'     the training set are left unchanged and a warning is issued.
#'     The training means and SDs are stored for application to the test set.
#' }
#'
#' **Step 2 - Geneset aggregation** collapses groups of features into single
#' summary features. If \code{genesets} is provided, features not present in
#' any geneset are discarded. Aggregation is performed on the (possibly
#' transformed) output of Step 1. Supported aggregation methods:
#' \itemize{
#'   \item \code{"mean"}: mean expression across geneset members.
#'   \item \code{"median"}: median expression across geneset members.
#'   \item \code{"sum"}: sum of expression across geneset members.
#'   \item \code{"pc1"}: first principal component of the geneset members,
#'     computed by PCA on the training samples. The feature loadings are stored
#'     and applied to the test set via \code{\link{predict_engineering}}.
#' }
#' If \code{genesets = NULL}, Step 2 is skipped and the output of Step 1 is
#' returned directly.
#'
#' @param X_train Numeric matrix of dimensions n (samples) x p (features).
#'   Training predictor matrix. Column names must be present and are used to
#'   match features to genesets. When \code{fc = TRUE}, must contain both
#'   pre-treatment and post-treatment rows for training individuals.
#' @param params A named list of engineering parameters with the following
#'   elements:
#'   \describe{
#'     \item{\code{method}}{Character string. Must be \code{"engineer"}.
#'       Required by the predictomics pipeline convention.}
#'     \item{\code{fc}}{Logical. If \code{TRUE}, fold-change engineering
#'       (Step 0) is applied before column transformation and geneset
#'       aggregation. Requires \code{treatment} and \code{participant_id}
#'       to be supplied to \code{run_engineering}. Defaults to \code{FALSE}.}
#'     \item{\code{col_transform}}{Character string. Column-wise transformation
#'       to apply. One of \code{"none"} (default) or \code{"z"} (z-score).}
#'     \item{\code{genesets}}{Named list of character vectors, or \code{NULL}
#'       (default). Each element is a geneset: a character vector of feature
#'       names corresponding to column names of \code{X_train}. Features not
#'       present in any geneset are discarded. Pass \code{NULL} to skip
#'       geneset aggregation.}
#'     \item{\code{agg_method}}{Character string. Aggregation method to apply
#'       within each geneset. One of \code{"mean"}, \code{"median"},
#'       \code{"sum"}, or \code{"pc1"}. Required when \code{genesets} is not
#'       \code{NULL}.}
#'   }
#' @param treatment Binary numeric vector (0/1) of length n. Required when
#'   \code{params$fc = TRUE}. Values of 0 indicate pre-treatment and 1
#'   indicate post-treatment. Pass \code{NULL} (default) otherwise.
#' @param participant_id Vector of length n identifying which rows belong to
#'   the same individual. Required when \code{params$fc = TRUE}. Each
#'   individual must have exactly one pre-treatment (0) and one post-treatment
#'   (1) row. Pass \code{NULL} (default) otherwise.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{X_transformed}}{Numeric matrix. The engineered training
#'       matrix. If FC engineering was applied, has \code{n_individuals} rows.
#'       If geneset aggregation was performed, columns correspond to genesets;
#'       otherwise columns correspond to the (transformed) input features.}
#'     \item{\code{fit}}{A named list of fitted parameters required to apply
#'       the same transformations to a test matrix via
#'       \code{\link{predict_engineering}}. Contains: \code{fc},
#'       \code{col_transform}, \code{col_means}, \code{col_sds} (for z-score),
#'       \code{genesets}, \code{agg_method}, and \code{pc1_loadings} (for
#'       pc1 aggregation).}
#'   }
#'
#' @seealso \code{\link{predict_engineering}}, \code{\link{predict_cv}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Simulate paired pre/post data: 20 individuals, 2 timepoints, 10 genes
#' n_ind <- 20
#' X <- matrix(rnorm(n_ind * 2 * 10), nrow = n_ind * 2, ncol = 10)
#' colnames(X) <- paste0("gene", 1:10)
#' treatment      <- rep(c(0, 1), each = n_ind)
#' participant_id <- rep(seq_len(n_ind), times = 2)
#'
#' # FC engineering only
#' params <- list(method = "engineer", fc = TRUE)
#' result <- run_engineering(X_train = X, params = params,
#'                           treatment = treatment,
#'                           participant_id = participant_id)
#' dim(result$X_transformed)  # 20 x 10
#'
#' # FC + z-score + geneset aggregation
#' genesets <- list(setA = paste0("gene", 1:5), setB = paste0("gene", 6:10))
#' params <- list(method = "engineer", fc = TRUE, col_transform = "z",
#'                genesets = genesets, agg_method = "mean")
#' result <- run_engineering(X_train = X, params = params,
#'                           treatment = treatment,
#'                           participant_id = participant_id)
#' dim(result$X_transformed)  # 20 x 2
#' }
#'
#' @export
# -----------------------------------------------------------------------------
run_engineering <- function(X_train, params, treatment = NULL,
                            participant_id = NULL) {

  # ---------------------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------------------
  .validate_X(X_train)
  .validate_engineering_params(params)

  if (is.null(colnames(X_train)))
    stop("[predictomics] X_train must have column names for feature engineering.",
         call. = FALSE)

  fc            <- isTRUE(params$fc)
  col_transform <- params$col_transform %||% "none"
  genesets      <- params$genesets
  agg_method    <- params$agg_method

  if (fc) {
    if (is.null(treatment) || is.null(participant_id))
      stop("[predictomics] treatment and participant_id must be supplied to ",
           "run_engineering() when params$fc = TRUE.", call. = FALSE)
    .validate_fc_inputs(treatment, participant_id, X_train)
  }

  # ---------------------------------------------------------------------------
  # 2. Step 0 - Fold-change engineering
  # ---------------------------------------------------------------------------
  if (fc) {
    X_train <- .compute_fc(X_train, treatment, participant_id)
  }

  # ---------------------------------------------------------------------------
  # 3. Step 1 - Column-wise transformation
  # ---------------------------------------------------------------------------
  col_means <- NULL
  col_sds   <- NULL

  X_out <- switch(col_transform,

                  none = X_train,

                  z = {
                    col_means <- colMeans(X_train)
                    col_sds   <- apply(X_train, 2, sd)

                    zero_var <- col_sds == 0
                    if (any(zero_var)) {
                      warning(
                        "[predictomics] ", sum(zero_var), " feature(s) have zero variance ",
                        "in the training set and will not be z-scored: ",
                        paste(colnames(X_train)[zero_var], collapse = ", "),
                        call. = FALSE
                      )
                      col_sds[zero_var] <- 1
                    }

                    sweep(sweep(X_train, 2, col_means, "-"), 2, col_sds, "/")
                  }
  )

  # ---------------------------------------------------------------------------
  # 4. Step 2 - Geneset aggregation
  # ---------------------------------------------------------------------------
  pc1_loadings <- NULL

  if (!is.null(genesets)) {

    feature_names <- colnames(X_out)

    X_out <- .aggregate_genesets(
      X             = X_out,
      genesets      = genesets,
      agg_method    = agg_method,
      feature_names = feature_names,
      is_train      = TRUE,
      pc1_loadings  = NULL,
      col_transform = col_transform
    )

    pc1_loadings <- attr(X_out, "pc1_loadings")
    attr(X_out, "pc1_loadings") <- NULL
  }

  # ---------------------------------------------------------------------------
  # 5. Assemble and return
  # ---------------------------------------------------------------------------
  fit <- list(
    fc            = fc,
    col_transform = col_transform,
    col_means     = col_means,
    col_sds       = col_sds,
    genesets      = genesets,
    agg_method    = agg_method,
    pc1_loadings  = pc1_loadings
  )

  list(X_transformed = X_out, fit = fit)
}


# -----------------------------------------------------------------------------
#' Apply fitted feature engineering to a new (test) matrix
#'
#' @description
#' Applies the feature engineering transformations fitted by
#' \code{\link{run_engineering}} on training data to a new matrix \code{X_new},
#' using only the parameters stored in \code{fit}. No parameters are
#' re-estimated from \code{X_new}.
#'
#' @param fit A fit object returned in the \code{fit} element of
#'   \code{\link{run_engineering}}.
#' @param X_new Numeric matrix. New predictor matrix to transform. Must have
#'   the same column names as the training matrix passed to
#'   \code{\link{run_engineering}}. When \code{fit$fc = TRUE}, must contain
#'   both pre-treatment and post-treatment rows for test individuals.
#' @param treatment Binary numeric vector (0/1). Required when
#'   \code{fit$fc = TRUE}. Same length as \code{nrow(X_new)}.
#' @param participant_id Vector identifying individual membership. Required
#'   when \code{fit$fc = TRUE}. Same length as \code{nrow(X_new)}.
#'
#' @return A numeric matrix of transformed features.
#'
#' @seealso \code{\link{run_engineering}}
#'
#' @export
# -----------------------------------------------------------------------------
predict_engineering <- function(fit, X_new, treatment = NULL,
                                participant_id = NULL) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!is.matrix(X_new) || !is.numeric(X_new))
    stop("[predictomics] X_new must be a numeric matrix.", call. = FALSE)
  if (is.null(colnames(X_new)))
    stop("[predictomics] X_new must have column names.", call. = FALSE)

  if (isTRUE(fit$fc)) {
    if (is.null(treatment) || is.null(participant_id))
      stop("[predictomics] treatment and participant_id must be supplied to ",
           "predict_engineering() when fit$fc = TRUE.", call. = FALSE)
    .validate_fc_inputs(treatment, participant_id, X_new)
  }

  # ---------------------------------------------------------------------------
  # 2. Step 0 - Apply fold-change engineering (parameter-free)
  # ---------------------------------------------------------------------------
  if (isTRUE(fit$fc)) {
    X_new <- .compute_fc(X_new, treatment, participant_id)
  }

  # ---------------------------------------------------------------------------
  # 3. Step 1 - Apply column-wise transformation using training parameters
  # ---------------------------------------------------------------------------
  X_out <- switch(fit$col_transform,

                  none = X_new,

                  z = sweep(sweep(X_new, 2, fit$col_means, "-"), 2, fit$col_sds, "/")
  )

  # ---------------------------------------------------------------------------
  # 4. Step 2 - Apply geneset aggregation using training parameters
  # ---------------------------------------------------------------------------
  if (!is.null(fit$genesets)) {
    X_out <- .aggregate_genesets(
      X             = X_out,
      genesets      = fit$genesets,
      agg_method    = fit$agg_method,
      feature_names = colnames(X_out),
      is_train      = FALSE,
      pc1_loadings  = fit$pc1_loadings,
      col_transform = fit$col_transform
    )
  }

  X_out
}


# =============================================================================
# Internal helpers
# =============================================================================

# -----------------------------------------------------------------------------
#' Compute fold-change matrix from paired pre/post data
#'
#' @description
#' For each individual and each feature, computes the difference
#' \eqn{X_{post} - X_{pre}}. Returns a matrix with one row per individual,
#' ordered by the sorted unique values of \code{participant_id}.
#'
#' @param X Numeric matrix with both pre (treatment == 0) and post
#'   (treatment == 1) rows.
#' @param treatment Binary numeric vector (0/1).
#' @param participant_id Vector of individual identifiers.
#'
#' @return Numeric matrix of dimensions n_individuals x p, with row names
#'   set to the sorted unique participant IDs.
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_fc <- function(X, treatment, participant_id) {

  ids      <- sort(unique(participant_id))
  n_ind    <- length(ids)
  p        <- ncol(X)
  X_fc     <- matrix(NA_real_, nrow = n_ind, ncol = p,
                     dimnames = list(as.character(ids), colnames(X)))

  for (i in seq_len(n_ind)) {
    id       <- ids[i]
    post_row <- which(participant_id == id & treatment == 1)
    pre_row  <- which(participant_id == id & treatment == 0)
    X_fc[i, ] <- X[post_row, ] - X[pre_row, ]
  }

  X_fc
}


# -----------------------------------------------------------------------------
#' Aggregate features into genesets
#'
#' @description
#' Internal workhorse for geneset aggregation. Called by both
#' \code{run_engineering} (training) and \code{predict_engineering} (test).
#' When \code{is_train = TRUE} and \code{agg_method = "pc1"}, PC1 loadings
#' are computed and attached as an attribute of the returned matrix for
#' retrieval by \code{run_engineering}.
#'
#' @param X Numeric matrix post column-wise transformation.
#' @param genesets Named list of character vectors of feature names.
#' @param agg_method Character string. One of "mean", "median", "sum", "pc1".
#' @param feature_names Character vector of column names of \code{X}.
#' @param is_train Logical. If \code{TRUE}, PC1 loadings are fitted from
#'   \code{X}. If \code{FALSE}, \code{pc1_loadings} must be supplied.
#' @param pc1_loadings Named list of PC1 loading vectors (one per geneset),
#'   or \code{NULL} when \code{is_train = TRUE}.
#' @param col_transform Character string. The column transformation applied
#'   prior to aggregation, used to determine centering/scaling in PCA.
#'
#' @return Numeric matrix of aggregated features (samples x genesets), with
#'   PC1 loadings attached as an attribute when \code{is_train = TRUE} and
#'   \code{agg_method = "pc1"}.
#'
#' @keywords internal
# -----------------------------------------------------------------------------
.aggregate_genesets <- function(X, genesets, agg_method, feature_names,
                                is_train, pc1_loadings, col_transform) {

  n_sets       <- length(genesets)
  n_samples    <- nrow(X)
  set_names    <- names(genesets)
  X_agg        <- matrix(NA_real_, nrow = n_samples, ncol = n_sets,
                         dimnames = list(rownames(X), set_names))
  pc1_loadings_out <- if (is_train && agg_method == "pc1")
    vector("list", n_sets) else NULL
  if (!is.null(pc1_loadings_out)) names(pc1_loadings_out) <- set_names

  for (i in seq_len(n_sets)) {

    gs_name  <- set_names[i]
    gs_genes <- intersect(genesets[[gs_name]], feature_names)
    X_sub    <- X[, gs_genes, drop = FALSE]

    if (agg_method == "pc1") {
      if (is_train) {
        do_scale <- col_transform != "z"
        pca      <- prcomp(X_sub, center = do_scale, scale. = do_scale)
        loadings <- pca$rotation[, 1]
        pc1_loadings_out[[gs_name]] <- list(
          loadings = loadings,
          center   = if (do_scale) pca$center else NULL,
          scale    = if (do_scale) pca$scale  else NULL
        )
        X_agg[, i] <- as.numeric(X_sub %*% loadings)
      } else {
        ls       <- pc1_loadings[[gs_name]]
        X_sub_sc <- if (!is.null(ls$center)) {
          scale(X_sub, center = ls$center, scale = ls$scale)
        } else {
          X_sub
        }
        X_agg[, i] <- as.numeric(X_sub_sc %*% ls$loadings)
      }
    } else {
      X_agg[, i] <- switch(agg_method,
                           mean   = rowMeans(X_sub),
                           median = apply(X_sub, 1, median),
                           sum    = rowSums(X_sub)
      )
    }
  }

  if (!is.null(pc1_loadings_out))
    attr(X_agg, "pc1_loadings") <- pc1_loadings_out

  X_agg
}

