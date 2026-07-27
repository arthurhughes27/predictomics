# =============================================================================
# run_engineering.R
# Feature engineering for the predictomics package.
#
# Implements two sequential steps:
#   1. Column-wise transformation of individual features (z-score or none).
#   2. Geneset aggregation: summarising groups of features into single values
#      (mean, median, sum, or pc1).
#
# Step 1 always precedes step 2. Both steps store the parameters fitted on
# the training set so that they can be applied to the test set without
# data leakage via predict_engineering().
# =============================================================================


# -----------------------------------------------------------------------------
#' Apply feature engineering to a training matrix
#'
#' @description
#' Applies a two-step feature engineering pipeline to the predictor matrix
#' \code{X_train}: (1) optional column-wise transformation of individual
#' features, followed by (2) optional aggregation of features into genesets.
#' All parameters required to apply the same transformations to a test set
#' are stored in the returned fit object and consumed by
#' \code{\link{predict_engineering}}.
#'
#' @details
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
#'   \item \code{"ssgsea"}: single-sample Gene Set Enrichment Analysis scores,
#'     computed via \code{GSVA::gsva()} with \code{GSVA::ssgseaParam()}
#'     (requires the Suggested \pkg{GSVA} package). Unlike the other
#'     aggregation methods, ssGSEA scores are recomputed independently on
#'     whatever matrix is passed to \code{run_engineering}/
#'     \code{predict_engineering} (there are no training-derived parameters to
#'     store, since the enrichment statistic is computed per sample). See
#'     \code{ssgsea_alpha}, \code{ssgsea_min_size}, \code{ssgsea_max_size}, and
#'     \code{ssgsea_normalize} below for tuning options. Because ssGSEA relies
#'     on the within-sample ranking of genes, \code{col_transform = "z"}
#'     (which rescales each gene/column differently) alters that ranking and
#'     is discouraged in combination with \code{agg_method = "ssgsea"}; a
#'     warning is issued if both are used together.
#' }
#' If \code{genesets = NULL}, Step 2 is skipped and the output of Step 1 is
#' returned directly.
#'
#' @param X_train Numeric matrix of dimensions n (samples) x p (features).
#'   Training predictor matrix. Column names must be present and are used to
#'   match features to genesets.
#' @param params A named list of engineering parameters with the following
#'   elements:
#'   \describe{
#'     \item{\code{method}}{Character string. Must be \code{"engineer"}.
#'       Required by the predictomics pipeline convention.}
#'     \item{\code{col_transform}}{Character string. Column-wise transformation
#'       to apply. One of \code{"none"} (default) or \code{"z"} (z-score).}
#'     \item{\code{genesets}}{Named list of character vectors, or \code{NULL}
#'       (default). Each element is a geneset: a character vector of feature
#'       names corresponding to column names of \code{X_train}. Features not
#'       present in any geneset are discarded. Pass \code{NULL} to skip
#'       geneset aggregation.}
#'     \item{\code{agg_method}}{Character string. Aggregation method to apply
#'       within each geneset. One of \code{"mean"}, \code{"median"},
#'       \code{"sum"}, \code{"pc1"}, or \code{"ssgsea"}. Required when
#'       \code{genesets} is not \code{NULL}.}
#'     \item{\code{ssgsea_alpha}}{Numeric. The exponent weight applied to the
#'       running-sum statistic in ssGSEA. Defaults to \code{0.25}. Only used
#'       when \code{agg_method = "ssgsea"}.}
#'     \item{\code{ssgsea_min_size}}{Positive integer. Minimum geneset size
#'       (after intersecting with the features present in \code{X_train})
#'       required for a geneset to be scored. Defaults to \code{1}. Only used
#'       when \code{agg_method = "ssgsea"}.}
#'     \item{\code{ssgsea_max_size}}{Numeric. Maximum geneset size allowed.
#'       Defaults to \code{Inf}. Only used when \code{agg_method = "ssgsea"}.}
#'     \item{\code{ssgsea_normalize}}{Logical. Whether to apply the Barbie et
#'       al. (2009) normalisation of ssGSEA scores. Defaults to \code{TRUE}.
#'       Only used when \code{agg_method = "ssgsea"}.}
#'   }
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{X_transformed}}{Numeric matrix. The engineered training
#'       matrix, with samples as rows. If geneset aggregation was performed,
#'       columns correspond to genesets; otherwise columns correspond to the
#'       (transformed) input features.}
#'     \item{\code{fit}}{A named list of fitted parameters required to apply
#'       the same transformations to a test matrix via
#'       \code{\link{predict_engineering}}. Contains:
#'       \code{col_transform}, \code{col_means}, \code{col_sds} (for z-score),
#'       \code{genesets}, \code{agg_method}, \code{pc1_loadings} (for pc1
#'       aggregation), and \code{ssgsea_alpha}, \code{ssgsea_min_size},
#'       \code{ssgsea_max_size}, \code{ssgsea_normalize} (for ssgsea
#'       aggregation).}
#'   }
#'
#' @seealso \code{\link{predict_engineering}}, \code{\link{predict_cv}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
#' colnames(X) <- paste0("gene", 1:20)
#'
#' genesets <- list(
#'   setA = paste0("gene", 1:10),
#'   setB = paste0("gene", 11:20)
#' )
#'
#' # Z-score genes then aggregate by geneset mean
#' params <- list(method = "engineer", col_transform = "z",
#'                genesets = genesets, agg_method = "mean")
#' result <- run_engineering(X_train = X, params = params)
#' dim(result$X_transformed)  # 50 x 2
#'
#' # Apply same transformations to test set
#' X_test <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
#' colnames(X_test) <- paste0("gene", 1:20)
#' X_test_transformed <- predict_engineering(result$fit, X_new = X_test)
#' }
#'
#' @export
# -----------------------------------------------------------------------------
run_engineering <- function(X_train, params) {

  # ---------------------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------------------
  .validate_X(X_train)
  .validate_engineering_params(params)

  if (is.null(colnames(X_train))) {
    stop("[predictomics] X_train must have column names for feature engineering.",
         call. = FALSE)
  }

  col_transform    <- params$col_transform %||% "none"
  genesets         <- params$genesets
  agg_method       <- params$agg_method
  ssgsea_alpha     <- params$ssgsea_alpha     %||% 0.25
  ssgsea_min_size  <- params$ssgsea_min_size  %||% 1L
  ssgsea_max_size  <- params$ssgsea_max_size  %||% Inf
  ssgsea_normalize <- params$ssgsea_normalize %||% TRUE

  if (identical(agg_method, "ssgsea") && identical(col_transform, "z"))
    warning(
      "[predictomics] col_transform = 'z' rescales each feature/column ",
      "differently, which alters each sample's within-sample gene ranking ",
      "and therefore the ssGSEA enrichment scores. Consider ",
      "col_transform = 'none' when agg_method = 'ssgsea'.",
      call. = FALSE
    )

  # ---------------------------------------------------------------------------
  # 2. Step 1 - Column-wise transformation
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
        col_sds[zero_var] <- 1  # avoid division by zero; feature stays as-is
      }

      sweep(sweep(X_train, 2, col_means, "-"), 2, col_sds, "/")
    }
  )

  # ---------------------------------------------------------------------------
  # 3. Step 2 - Geneset aggregation
  # ---------------------------------------------------------------------------
  pc1_loadings <- NULL

  if (!is.null(genesets)) {

    feature_names <- colnames(X_out)

    X_out <- .aggregate_genesets(
      X                = X_out,
      genesets         = genesets,
      agg_method       = agg_method,
      feature_names    = feature_names,
      is_train         = TRUE,
      pc1_loadings     = NULL,
      col_transform    = col_transform,
      ssgsea_alpha     = ssgsea_alpha,
      ssgsea_min_size  = ssgsea_min_size,
      ssgsea_max_size  = ssgsea_max_size,
      ssgsea_normalize = ssgsea_normalize
    )

    pc1_loadings <- attr(X_out, "pc1_loadings")
    attr(X_out, "pc1_loadings") <- NULL
  }

  # ---------------------------------------------------------------------------
  # 4. Assemble and return
  # ---------------------------------------------------------------------------
  fit <- list(
    col_transform    = col_transform,
    col_means        = col_means,
    col_sds          = col_sds,
    genesets         = genesets,
    agg_method       = agg_method,
    pc1_loadings     = pc1_loadings,
    ssgsea_alpha     = ssgsea_alpha,
    ssgsea_min_size  = ssgsea_min_size,
    ssgsea_max_size  = ssgsea_max_size,
    ssgsea_normalize = ssgsea_normalize
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
#'   \code{\link{run_engineering}}.
#'
#' @return A numeric matrix of transformed features, with the same number of
#'   rows as \code{X_new}.
#'
#' @seealso \code{\link{run_engineering}}
#'
#' @export
# -----------------------------------------------------------------------------
predict_engineering <- function(fit, X_new) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!is.matrix(X_new) || !is.numeric(X_new))
    stop("[predictomics] X_new must be a numeric matrix.", call. = FALSE)
  if (is.null(colnames(X_new)))
    stop("[predictomics] X_new must have column names.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2. Step 1 - Apply column-wise transformation using training parameters
  # ---------------------------------------------------------------------------
  X_out <- switch(fit$col_transform,

    none = X_new,

    z = sweep(sweep(X_new, 2, fit$col_means, "-"), 2, fit$col_sds, "/")
  )

  # ---------------------------------------------------------------------------
  # 3. Step 2 - Apply geneset aggregation using training parameters
  # ---------------------------------------------------------------------------
  if (!is.null(fit$genesets)) {
    X_out <- .aggregate_genesets(
      X                = X_out,
      genesets         = fit$genesets,
      agg_method       = fit$agg_method,
      feature_names    = colnames(X_out),
      is_train         = FALSE,
      pc1_loadings     = fit$pc1_loadings,
      col_transform    = fit$col_transform,
      ssgsea_alpha     = fit$ssgsea_alpha,
      ssgsea_min_size  = fit$ssgsea_min_size,
      ssgsea_max_size  = fit$ssgsea_max_size,
      ssgsea_normalize = fit$ssgsea_normalize
    )
  }

  X_out
}


# =============================================================================
# Internal helpers
# =============================================================================

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
#' @param agg_method Character string. One of "mean", "median", "sum", "pc1",
#'   "ssgsea".
#' @param feature_names Character vector of column names of \code{X}.
#' @param is_train Logical. If \code{TRUE}, PC1 loadings are fitted from
#'   \code{X}. If \code{FALSE}, \code{pc1_loadings} must be supplied. Ignored
#'   for \code{agg_method = "ssgsea"}, which is recomputed from \code{X}
#'   regardless.
#' @param pc1_loadings Named list of PC1 loading vectors (one per geneset),
#'   or \code{NULL} when \code{is_train = TRUE}.
#' @param ssgsea_alpha,ssgsea_min_size,ssgsea_max_size,ssgsea_normalize
#'   ssGSEA tuning parameters, only used when \code{agg_method = "ssgsea"}.
#'   See \code{\link{run_engineering}}.
#'
#' @return Numeric matrix of aggregated features (samples x genesets), with
#'   PC1 loadings attached as an attribute when \code{is_train = TRUE} and
#'   \code{agg_method = "pc1"}.
#'
#' @keywords internal
# -----------------------------------------------------------------------------
.aggregate_genesets <- function(X, genesets, agg_method, feature_names,
                                is_train, pc1_loadings, col_transform,
                                ssgsea_alpha = 0.25, ssgsea_min_size = 1L,
                                ssgsea_max_size = Inf, ssgsea_normalize = TRUE) {

  if (agg_method == "ssgsea") {
    return(.compute_ssgsea(
      X         = X,
      genesets  = genesets,
      alpha     = ssgsea_alpha,
      min_size  = ssgsea_min_size,
      max_size  = ssgsea_max_size,
      normalize = ssgsea_normalize
    ))
  }

  n_sets       <- length(genesets)
  n_samples    <- nrow(X)
  set_names    <- names(genesets)
  X_agg        <- matrix(NA_real_, nrow = n_samples, ncol = n_sets,
                          dimnames = list(rownames(X), set_names))
  pc1_loadings_out <- if (is_train && agg_method == "pc1") vector("list", n_sets) else NULL
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


# -----------------------------------------------------------------------------
#' Compute single-sample GSEA (ssGSEA) enrichment scores
#'
#' @description
#' Internal workhorse for \code{agg_method = "ssgsea"}, called by
#' \code{\link{.aggregate_genesets}}. Computes ssGSEA scores via
#' \code{GSVA::gsva()}/\code{GSVA::ssgseaParam()}. Unlike the other
#' aggregation methods, ssGSEA scores are computed independently on whatever
#' matrix is supplied (there are no training-derived parameters that can be
#' stored and reapplied, since the enrichment statistic is computed per
#' sample from that sample's own gene ranking).
#'
#' @param X Numeric matrix of dimensions n (samples) x p (features), post
#'   column-wise transformation.
#' @param genesets Named list of character vectors of feature names.
#' @param alpha Numeric. Exponent weight for the running-sum statistic.
#' @param min_size Positive integer. Minimum geneset size (after intersecting
#'   with \code{colnames(X)}) required for a geneset to be scored.
#' @param max_size Numeric. Maximum geneset size allowed.
#' @param normalize Logical. Whether to apply the Barbie et al. (2009)
#'   normalisation of ssGSEA scores.
#'
#' @return Numeric matrix of ssGSEA scores (samples x genesets), in the same
#'   geneset order as \code{genesets}.
#'
#' @keywords internal
# -----------------------------------------------------------------------------
.compute_ssgsea <- function(X, genesets, alpha, min_size, max_size, normalize) {

  if (!requireNamespace("GSVA", quietly = TRUE))
    stop(
      "[predictomics] Package 'GSVA' is required for agg_method = 'ssgsea'. ",
      "Please install it (e.g. via BiocManager::install('GSVA')).",
      call. = FALSE
    )

  feature_names     <- colnames(X)
  genesets_filtered <- lapply(genesets, intersect, y = feature_names)

  empty <- vapply(genesets_filtered, length, integer(1)) == 0L
  if (any(empty))
    stop(
      "[predictomics] The following geneset(s) have no overlapping features ",
      "with X and cannot be scored via ssgsea: ",
      paste(names(genesets_filtered)[empty], collapse = ", "), ".",
      call. = FALSE
    )

  expr <- t(X)  # genes (features) x samples

  par <- GSVA::ssgseaParam(
    exprData  = expr,
    geneSets  = genesets_filtered,
    minSize   = min_size,
    maxSize   = max_size,
    alpha     = alpha,
    normalize = normalize
  )

  scores <- GSVA::gsva(par)
  X_agg  <- t(scores)
  X_agg  <- X_agg[, names(genesets_filtered), drop = FALSE]
  rownames(X_agg) <- rownames(X)

  X_agg
}
