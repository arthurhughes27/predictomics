# =============================================================================
# print.R
# Print method for predictomics result objects.
# =============================================================================


# -----------------------------------------------------------------------------
#' Print a predictomics result object
#'
#' @description
#' Prints a concise, human-readable summary of a \code{predictomics} result
#' object to the console, including the call, dataset dimensions, pipeline
#' configuration, and performance metrics.
#'
#' @param x A \code{predictomics} object returned by \code{\link{predict_cv}}.
#' @param digits Integer. Number of decimal places for metrics. Defaults to
#'   \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{predict_cv}}, \code{\link{metrics.predictomics}},
#'   \code{\link{plot.predictomics}}
#'
#' @export
# -----------------------------------------------------------------------------
print.predictomics <- function(x, digits = 4, ...) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!inherits(x, "predictomics"))
    stop("[predictomics] x must be a predictomics object returned by ",
         "predict_cv().", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2. Resolve pipeline labels
  # ---------------------------------------------------------------------------
  eng_label <- if (is.null(x$engineering_params)) {
    "none"
  } else {
    paste0(x$engineering_params$method,
           if (!is.null(x$engineering_params$col_transform))
             paste0(" [col_transform = ", x$engineering_params$col_transform, "]"),
           if (!is.null(x$engineering_params$genesets))
             paste0(" [genesets: ", length(x$engineering_params$genesets),
                    " sets, agg = ", x$engineering_params$agg_method, "]"))
  }

  sel_label <- if (is.null(x$selection_params)) {
    "none"
  } else {
    x$selection_params$method
  }

  mod_label <- x$model_params$method

  cv_label  <- if (!is.null(x$cv_type) && x$cv_type == "loo") {
    "leave-one-out"
  } else {
    paste0(x$n_folds %||% "?", "-fold")
  }

  # ---------------------------------------------------------------------------
  # 3. Compute metrics
  # ---------------------------------------------------------------------------
  m <- metrics.predictomics(x, digits = digits)

  # ---------------------------------------------------------------------------
  # 4. Print
  # ---------------------------------------------------------------------------
  cat("\n")
  cat("=================================================\n")
  cat(" predictomics: cross-validated prediction result \n")
  cat("=================================================\n")

  cat("\nCall:\n ")
  cat(deparse(x$call), "\n")

  cat("\nData:\n")
  cat("  Samples  :", x$n_samples, "\n")
  cat("  Features :", x$n_features_input, "\n")

  cat("\nCross-validation:\n")
  cat("  Type     :", cv_label, "\n")
  if (!is.null(x$outside_cv) && x$outside_cv)
    cat("  WARNING  : outside_cv = TRUE (estimates may be biased)\n")

  cat("\nPipeline:\n")
  cat("  Engineering :", eng_label, "\n")
  cat("  Selection   :", sel_label, "\n")
  cat("  Model       :", mod_label, "\n")

  cat("\nPerformance metrics:\n")
  cat("  RMSE      :", m["RMSE"],      "\n")
  cat("  sRMSE     :", m["sRMSE"],     "\n")
  cat("  R2        :", m["R2"],        "\n")
  cat("  SpearmanR :", m["SpearmanR"], "\n")

  cat("=================================================\n\n")

  invisible(x)
}
