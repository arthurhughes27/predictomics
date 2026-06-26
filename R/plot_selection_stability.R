# =============================================================================
# plot_selection_stability.R
# Visualisation of variable selection stability across CV folds for
# predictomics result objects.
# =============================================================================


# -----------------------------------------------------------------------------
#' Plot variable selection stability across cross-validation folds
#'
#' @description
#' Visualises the stability of variable selection across CV folds from a
#' \code{predictomics} result object. Two complementary plots are available:
#' a selection frequency bar chart and a binary selection heatmap. Both are
#' based on the per-fold selection diagnostics stored in the result object.
#'
#' @details
#' **Selection frequency bar chart**: for each feature selected in at least
#' one fold, shows the proportion of folds in which it was selected. Features
#' are sorted by selection frequency (descending). A dashed reference line at
#' \eqn{1/k} marks the random selection frequency; a second line at 1.0 marks
#' features selected in every fold.
#'
#' **Selection heatmap**: a binary grid with features on the y-axis (sorted
#' by selection frequency, most stable at top) and fold indices on the x-axis.
#' A filled cell indicates the feature was selected in that fold. This reveals
#' whether instability is driven by specific folds or is diffuse across all
#' folds.
#'
#' Both plots show the top \code{top_n} features by selection frequency. For
#' LOO CV with many folds, the heatmap is automatically suppressed and only
#' the frequency bar chart is returned, with a message explaining this.
#'
#' @param x A \code{predictomics} object returned by \code{\link{predict_cv}}.
#' @param type Character string. Which selection diagnostics to use. One of
#'   \code{"explicit"} (from \code{fold_selection_diagnostics}, i.e. the
#'   filter or wrapper selection method) or \code{"embedded"} (from
#'   \code{fold_embedded_selection_diagnostics}, i.e. lasso/glmnet non-zero
#'   coefficients). Defaults to \code{"explicit"}.
#' @param plot_type Character string. Which plot(s) to produce. One of
#'   \code{"frequency"}, \code{"heatmap"}, or \code{"both"} (default).
#' @param top_n Positive integer. Number of top features (by selection
#'   frequency) to display. Defaults to \code{20}.
#' @param heatmap_fold_threshold Integer. Maximum number of folds for which
#'   the heatmap is shown. For LOO CV or any setting with more folds than this
#'   threshold, the heatmap is suppressed. Defaults to \code{30}.
#' @param ... Additional arguments passed to \code{ggplot2::theme}.
#'
#' @return If \code{plot_type = "both"}, a named list with elements
#'   \code{"frequency"} and \code{"heatmap"} (the heatmap may be \code{NULL}
#'   if suppressed). If \code{plot_type = "frequency"} or
#'   \code{plot_type = "heatmap"}, a single \code{ggplot} object.
#'
#' @seealso \code{\link{predict_cv}}, \code{\link{plot.predictomics}}
#'
#' @examples
#' \dontrun{
#' result <- predict_cv(
#'   Y = Y, X = X,
#'   selection_params = list(method = "pearson", top_n = 50),
#'   model_params     = list(method = "lm")
#' )
#'
#' # Both plots (default)
#' plots <- plot_selection_stability(result)
#' plots$frequency
#' plots$heatmap
#'
#' # Frequency bar chart only
#' plot_selection_stability(result, plot_type = "frequency")
#'
#' # Embedded selection stability (lasso)
#' result2 <- predict_cv(Y = Y, X = X,
#'                       model_params = list(method = "lasso"))
#' plot_selection_stability(result2, type = "embedded")
#' }
#'
#' @export
# -----------------------------------------------------------------------------
plot_selection_stability <- function(x,
                                     type                  = "explicit",
                                     plot_type             = "both",
                                     top_n                 = 20L,
                                     heatmap_fold_threshold = 30L,
                                     ...) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!inherits(x, "predictomics"))
    stop("[predictomics] x must be a predictomics object returned by ",
         "predict_cv().", call. = FALSE)

  if (!type %in% c("explicit", "embedded"))
    stop("[predictomics] type must be 'explicit' or 'embedded'.", call. = FALSE)

  if (!plot_type %in% c("frequency", "heatmap", "both"))
    stop("[predictomics] plot_type must be 'frequency', 'heatmap', or 'both'.",
         call. = FALSE)

  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1L)
    stop("[predictomics] top_n must be a positive integer.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2. Extract diagnostics
  # ---------------------------------------------------------------------------
  if (type == "explicit") {

    if (!is.null(x$outside_cv) && x$outside_cv)
      stop("[predictomics] Selection stability plots require inside-CV ",
           "selection (outside_cv = FALSE). When outside_cv = TRUE, selection ",
           "is performed once on the full dataset and stability cannot be ",
           "assessed.", call. = FALSE)

    diag_list <- x$fold_selection_diagnostics

    if (is.null(diag_list))
      stop("[predictomics] No explicit selection diagnostics found. Run ",
           "predict_cv() with a selection_params argument and outside_cv = ",
           "FALSE.", call. = FALSE)

    type_label <- paste0("Explicit selection (", x$selection_params$method, ")")

  } else {

    diag_list <- x$fold_embedded_selection_diagnostics

    if (is.null(diag_list))
      stop("[predictomics] No embedded selection diagnostics found. Run ",
           "predict_cv() with model_params$method set to 'lasso' or 'glmnet'.",
           call. = FALSE)

    type_label <- paste0("Embedded selection (", x$model_params$method, ")")
  }

  # ---------------------------------------------------------------------------
  # 3. Compute selection frequency across folds
  # ---------------------------------------------------------------------------
  n_folds      <- length(diag_list)
  all_features <- unique(unlist(lapply(diag_list, `[[`, "selected_features")))

  if (length(all_features) == 0L)
    stop("[predictomics] No features were selected in any fold.", call. = FALSE)

  # Binary selection matrix: features x folds
  sel_matrix <- matrix(
    0L,
    nrow     = length(all_features),
    ncol     = n_folds,
    dimnames = list(all_features, paste0("Fold ", seq_len(n_folds)))
  )

  for (k in seq_len(n_folds)) {
    sel_k <- diag_list[[k]]$selected_features
    if (!is.null(sel_k) && length(sel_k) > 0L)
      sel_matrix[sel_k, k] <- 1L
  }

  # Selection frequency per feature
  freq           <- rowMeans(sel_matrix)
  freq_sorted    <- sort(freq, decreasing = TRUE)
  top_features   <- names(freq_sorted)[seq_len(min(top_n, length(freq_sorted)))]

  freq_df <- data.frame(
    feature   = factor(top_features, levels = rev(top_features)),
    frequency = freq_sorted[top_features]
  )

  # ---------------------------------------------------------------------------
  # 4. Frequency bar chart
  # ---------------------------------------------------------------------------
  p_freq <- ggplot2::ggplot(
    freq_df,
    ggplot2::aes(x = frequency, y = feature)
  ) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +

    # Reference line: random selection frequency
    # ggplot2::geom_vline(
    #   xintercept = 1 / n_folds,
    #   linetype   = "dashed",
    #   colour     = "grey50",
    #   linewidth  = 0.6
    # ) +
    # ggplot2::annotate(
    #   "text",
    #   x      = 1 / n_folds + 0.01,
    #   y      = 0.5,
    #   label  = paste0("1/k = ", round(1 / n_folds, 2)),
    #   hjust  = 0,
    #   size   = 3,
    #   colour = "grey40"
    # ) +

    # Reference line: selected in all folds
    # ggplot2::geom_vline(
    #   xintercept = 1,
    #   linetype   = "dotted",
    #   colour     = "darkgreen",
    #   linewidth  = 0.6
    # ) +

    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +

    ggplot2::labs(
      x        = "Selection frequency",
      y        = NULL,
      title    = "Variable selection stability",
      subtitle = paste0(type_label, "  |  Top ", nrow(freq_df),
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

  # ---------------------------------------------------------------------------
  # 5. Heatmap (suppressed for LOO or many-fold settings)
  # ---------------------------------------------------------------------------
  p_heat <- NULL

  if (plot_type %in% c("heatmap", "both")) {

    if (n_folds > heatmap_fold_threshold) {
      message(
        "[predictomics] Heatmap suppressed: number of folds (", n_folds,
        ") exceeds heatmap_fold_threshold (", heatmap_fold_threshold, "). ",
        "Use plot_type = 'frequency' or increase heatmap_fold_threshold."
      )
    } else {

      # Subset to top_n features, preserving frequency order
      heat_matrix <- sel_matrix[top_features, , drop = FALSE]

      heat_df <- as.data.frame(heat_matrix)
      heat_df$feature <- factor(rownames(heat_matrix),
                                levels = rev(top_features))
      heat_long <- tidyr::pivot_longer(
        heat_df,
        cols      = -feature,
        names_to  = "fold",
        values_to = "selected"
      )
      heat_long$fold     <- factor(heat_long$fold,
                                   levels = paste0("Fold ", seq_len(n_folds)))
      heat_long$selected <- factor(heat_long$selected,
                                   levels = c(0L, 1L),
                                   labels = c("Not selected", "Selected"))

      p_heat <- ggplot2::ggplot(
        heat_long,
        ggplot2::aes(x = fold, y = feature, fill = selected)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +

        ggplot2::scale_fill_manual(
          values = c("Not selected" = "grey92", "Selected" = "steelblue"),
          name   = NULL
        ) +

        ggplot2::labs(
          x        = "Fold",
          y        = NULL,
          title    = "Variable selection stability - fold detail",
          subtitle = paste0(type_label, "  |  Top ", nrow(freq_df),
                            " features  |  ", n_folds, " folds")
        ) +

        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title       = ggplot2::element_text(face = "bold", size = 12),
          plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 10),
          axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                                    size = 8),
          axis.text.y      = ggplot2::element_text(size = 9),
          panel.grid       = ggplot2::element_blank(),
          legend.position  = "right",
          ...
        )
    }
  }

  # ---------------------------------------------------------------------------
  # 6. Return
  # ---------------------------------------------------------------------------
  if (plot_type == "frequency") return(p_freq)
  if (plot_type == "heatmap")   return(p_heat)

  list(frequency = p_freq, heatmap = p_heat)
}
