# =============================================================================
# plot_feature_importance.R
# Visualisation of model-derived feature importance across CV folds for
# predictomics result objects.
# =============================================================================


# -----------------------------------------------------------------------------
#' Plot model-derived feature importance across cross-validation folds
#'
#' @description
#' Visualises per-feature importance, averaged across CV folds, from a
#' \code{predictomics} result object fitted with
#' \code{model_params$compute_importance = TRUE}. Unlike
#' \code{\link{plot_selection_stability}} (which shows how often a feature
#' was \emph{selected}), this shows the magnitude of each feature's
#' contribution to the fitted model, so it is meaningful even when no
#' explicit or embedded feature selection was used.
#'
#' @details
#' For each feature that reached the model in at least one fold, the mean
#' absolute importance across the folds in which it was present is computed
#' and features are ranked in decreasing order. The top \code{top_n} features
#' are shown as a horizontal bar chart, in the same visual style as
#' \code{\link{plot_selection_stability}}.
#'
#' Importance scores are method-specific (see \code{\link{run_model}}):
#' raw coefficients (signed) for \code{"lm"}, \code{"glmnet"}, \code{"lasso"},
#' \code{"ridge"}, and \code{"svr"}; impurity-based variable importance
#' (unsigned) for \code{"ranger"}. Because coefficient-based scores are only
#' comparable across features on a common scale, a warning is issued (in
#' addition to the one already issued by \code{\link{predict_cv}} at fit
#' time) when \code{x$model_params$scale} is not \code{TRUE} and the fitted
#' method is coefficient-based.
#'
#' @param x A \code{predictomics} object returned by \code{\link{predict_cv}},
#'   fitted with \code{model_params$compute_importance = TRUE}.
#' @param top_n Positive integer. Number of top features (by mean absolute
#'   importance) to display. Defaults to \code{20}.
#' @param ... Additional arguments passed to \code{ggplot2::theme}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{predict_cv}}, \code{\link{run_model}},
#'   \code{\link{plot_selection_stability}}
#'
#' @examples
#' \dontrun{
#' result <- predict_cv(
#'   Y = Y, X = X,
#'   model_params = list(method = "svr", scale = TRUE,
#'                       compute_importance = TRUE)
#' )
#' plot_feature_importance(result)
#' }
#'
#' @export
# -----------------------------------------------------------------------------
plot_feature_importance <- function(x, top_n = 20L, ...) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!inherits(x, "predictomics"))
    stop("[predictomics] x must be a predictomics object returned by ",
         "predict_cv().", call. = FALSE)

  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1L)
    stop("[predictomics] top_n must be a positive integer.", call. = FALSE)

  diag_list <- x$fold_feature_importance

  if (is.null(diag_list))
    stop("[predictomics] No feature importance diagnostics found. Run ",
         "predict_cv() with model_params$compute_importance = TRUE.",
         call. = FALSE)

  non_null_idx <- which(!vapply(diag_list, is.null, logical(1)))

  if (length(non_null_idx) == 0L)
    stop("[predictomics] No feature importance scores were recorded in any ",
         "fold.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2. Scale warning
  # ---------------------------------------------------------------------------
  importance_type <- diag_list[[non_null_idx[1]]]$type

  coefficient_based_methods <- c("lm", "glmnet", "ridge", "lasso", "svr")
  if (!is.null(importance_type) && importance_type == "coefficient" &&
      x$model_params$method %in% coefficient_based_methods &&
      !isTRUE(x$model_params$scale)) {
    warning(
      "[predictomics] Feature importance for method = '", x$model_params$method,
      "' is a raw (signed) coefficient, which is only comparable across ",
      "features that are on a common scale. This result object was fitted ",
      "with model_params$scale not TRUE, so the ranking below may be ",
      "misleading - refit with model_params$scale = TRUE for a meaningful ",
      "cross-feature comparison.",
      call. = FALSE, immediate. = TRUE
    )
  }

  # ---------------------------------------------------------------------------
  # 3. Compute mean absolute importance across folds
  # ---------------------------------------------------------------------------
  n_folds      <- length(diag_list)
  all_features <- unique(unlist(lapply(
    diag_list,
    function(d) { nm <- names(d$scores); nm[!is.na(nm)] }
  )))

  if (length(all_features) == 0L)
    stop("[predictomics] No feature importance scores were recorded in any ",
         "fold.", call. = FALSE)

  # features x folds, NA where the feature did not reach the model that fold
  imp_matrix <- matrix(
    NA_real_,
    nrow     = length(all_features),
    ncol     = n_folds,
    dimnames = list(all_features, paste0("Fold ", seq_len(n_folds)))
  )

  for (k in seq_len(n_folds)) {
    scores_k <- diag_list[[k]]$scores
    if (!is.null(scores_k) && length(scores_k) > 0L) {
      # Guard against any stray NA or unrecognised names rather than letting
      # them reach a matrix subscript (NA/unmatched names are never valid
      # row indices and would otherwise error out here).
      valid <- !is.na(names(scores_k)) & names(scores_k) %in% all_features
      if (any(valid))
        imp_matrix[names(scores_k)[valid], k] <- scores_k[valid]
    }
  }

  mean_abs_imp <- rowMeans(abs(imp_matrix), na.rm = TRUE)
  imp_sorted   <- sort(mean_abs_imp, decreasing = TRUE)
  top_features <- names(imp_sorted)[seq_len(min(top_n, length(imp_sorted)))]

  imp_df <- data.frame(
    feature    = factor(top_features, levels = rev(top_features)),
    importance = imp_sorted[top_features]
  )

  # ---------------------------------------------------------------------------
  # 4. Bar chart
  # ---------------------------------------------------------------------------
  type_label <- if (identical(importance_type, "impurity")) {
    "Impurity-based importance"
  } else {
    "Mean |coefficient|"
  }

  ggplot2::ggplot(
    imp_df,
    ggplot2::aes(x = importance, y = feature)
  ) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +

    ggplot2::labs(
      x        = type_label,
      y        = NULL,
      title    = "Feature importance",
      subtitle = paste0(x$model_params$method, "  |  Top ", nrow(imp_df),
                        " features  |  ", n_folds, " folds")
    ) +

    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 10),
      axis.text.y      = ggplot2::element_text(size = 9),
      axis.text.x      = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      ...
    )
}
