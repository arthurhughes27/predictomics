# =============================================================================
# plot.R
# Scatter plot of observed vs cross-validated predicted values for
# predictomics result objects.
# =============================================================================


# -----------------------------------------------------------------------------
#' Plot observed versus cross-validated predicted values
#'
#' @description
#' Produces a scatter plot of observed (x-axis) versus cross-validated
#' predicted (y-axis) values from a \code{predictomics} result object.
#' A diagonal reference line (y = x) representing perfect prediction is
#' included. sRMSE and R² are annotated in the top left corner of the plot.
#' If a \code{treatment} variable is present in the result object, points are
#' coloured by treatment group, unless \code{highlight} is supplied, in which
#' case it takes precedence.
#'
#' @details
#' Axis limits are computed as the range of the combined observed and predicted
#' values, expanded by 5\% on each side, and are identical for both axes so
#' that the aspect ratio is 1:1. This ensures that the y = x reference line
#' bisects the plot at 45 degrees and that deviations from perfect prediction
#' are visually unambiguous.
#'
#' R² is the squared Pearson correlation. sRMSE is RMSE divided by the
#' standard deviation of the observed values. Both are computed via
#' \code{\link{metrics.predictomics}}.
#'
#' When \code{highlight} is supplied, points are coloured by its values
#' instead of \code{treatment} (if any), and any treatment colouring is
#' overridden entirely. For \code{highlight_type = "categorical"}, a discrete
#' colour scale adapted to the number of levels is used - \code{RColorBrewer}
#' \code{"Set2"} for up to 8 levels, or \code{viridis} (discrete) beyond that
#' - with a legend showing each level. For \code{highlight_type =
#' "continuous"}, a \code{viridis} colour gradient is used, with a continuous
#' colourbar legend.
#'
#' When \code{highlight} is absent and \code{treatment} is present, a
#' colorblind-friendly palette from \code{RColorBrewer} (\code{"Set2"}) is
#' used to distinguish treatment groups. When neither is present, all points
#' are drawn in \code{point_colour}.
#'
#' @param x A \code{predictomics} object returned by \code{\link{predict_cv}}.
#' @param point_colour Character string. Colour of points when neither
#'   \code{highlight} nor a treatment variable is present. Defaults to
#'   \code{"steelblue"}.
#' @param point_alpha Numeric between 0 and 1. Transparency of points.
#'   Defaults to \code{0.6}.
#' @param point_size Numeric. Size of points. Defaults to \code{2}.
#' @param annotation_pos Character string. Position of the metric annotation.
#'   One of \code{"topleft"} (default) or \code{"bottomright"}.
#' @param highlight A vector used to colour points, overriding any treatment
#'   colouring. Must have the same length as \code{x$observed} (i.e. one
#'   value per plotted sample, in the same order). Pass \code{NULL} (default)
#'   to colour by treatment (if present) or not at all.
#' @param highlight_type Character string. Required when \code{highlight} is
#'   supplied. One of \code{"categorical"} (discrete colour scale, with a
#'   legend entry per level) or \code{"continuous"} (colour gradient with a
#'   colourbar legend).
#' @param highlight_label Character string. Legend title used for
#'   \code{highlight}. Defaults to \code{"Highlight"} if \code{NULL}.
#' @param ... Additional arguments passed to \code{ggplot2::theme}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{predict_cv}}, \code{\link{metrics.predictomics}},
#'   \code{\link{print.predictomics}}
#'
#' @examples
#' \dontrun{
#' result <- predict_cv(Y = Y, X = X, model_params = list(method = "glmnet"))
#' plot(result)
#'
#' # With treatment colouring
#' result <- predict_cv(Y = Y, X = X, treatment = treatment)
#' plot(result)
#'
#' # Highlight by a categorical variable (overrides treatment colouring)
#' plot(result, highlight = batch, highlight_type = "categorical",
#'     highlight_label = "Batch")
#'
#' # Highlight by a continuous variable
#' plot(result, highlight = age, highlight_type = "continuous",
#'     highlight_label = "Age")
#'
#' # Customise appearance
#' plot(result, point_colour = "darkred", point_alpha = 0.8)
#' }
#'
#' @export
# -----------------------------------------------------------------------------
plot.predictomics <- function(x,
                              point_colour    = "steelblue",
                              point_alpha     = 0.6,
                              point_size      = 2,
                              annotation_pos  = "topleft",
                              highlight       = NULL,
                              highlight_type  = NULL,
                              highlight_label = NULL,
                              ...) {

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!inherits(x, "predictomics"))
    stop("[predictomics] x must be a predictomics object returned by ",
         "predict_cv().", call. = FALSE)

  if (!annotation_pos %in% c("topleft", "bottomright"))
    stop("[predictomics] annotation_pos must be 'topleft' or 'bottomright'.",
         call. = FALSE)

  has_highlight <- !is.null(highlight)

  if (has_highlight) {

    if (length(highlight) != length(x$observed))
      stop("[predictomics] highlight must have length ",
           length(x$observed), " (one value per plotted sample, matching ",
           "x$observed).", call. = FALSE)

    if (is.null(highlight_type) ||
        !highlight_type %in% c("categorical", "continuous"))
      stop("[predictomics] highlight_type must be 'categorical' or ",
           "'continuous' when highlight is supplied.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 2. Compute metrics and axis limits
  # ---------------------------------------------------------------------------
  m        <- metrics.predictomics(x, digits = 3)
  obs      <- x$observed
  pred     <- x$predicted

  all_vals <- c(obs, pred)
  rng      <- range(all_vals, na.rm = TRUE)
  pad      <- diff(rng) * 0.05
  lims     <- c(rng[1] - pad, rng[2] + pad)

  # ---------------------------------------------------------------------------
  # 3. Annotation position
  # ---------------------------------------------------------------------------
  ann_x <- if (annotation_pos == "topleft") lims[1] + diff(lims) * 0.02
  else                             lims[2] - diff(lims) * 0.40
  ann_y <- if (annotation_pos == "topleft") lims[2] - diff(lims) * 0.02
  else                             lims[1] + diff(lims) * 0.12
  ann_vjust <- if (annotation_pos == "topleft") 1 else 0

  ann_label <- paste0(
    "sRMSE = ", m["sRMSE"], "\n",
    "R\u00B2 = ",  m["R2"]
  )

  # ---------------------------------------------------------------------------
  # 4. Prepare plot data frame and point colouring (highlight overrides
  # treatment entirely when supplied)
  # ---------------------------------------------------------------------------
  has_treatment <- !has_highlight && !is.null(x$treatment)

  plot_df <- data.frame(observed = obs, predicted = pred)

  if (has_highlight) {
    if (highlight_type == "categorical") {
      plot_df$highlight_value <- factor(highlight)
    } else {
      plot_df$highlight_value <- as.numeric(highlight)
    }
  } else if (has_treatment) {
    # Convert binary numeric to a labelled factor for the legend
    trt <- if (isTRUE(x$paired_rise)) {
      x$treatment[x$treatment == 1]
    } else {
      x$treatment
    }
    if (is.numeric(trt)) {
      trt <- factor(trt, levels = c(0, 1),
                    labels = c("Control (0)", "Active (1)"))
    }
    plot_df$treatment_group <- trt
  }

  use_colour  <- has_highlight || has_treatment
  colour_name <- highlight_label %||% "Highlight"

  # ---------------------------------------------------------------------------
  # 5. Build plot
  # ---------------------------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = observed, y = predicted)
  ) +

    # Reference line (perfect prediction)
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey40",
                         linewidth = 0.6) +

    # Points - coloured by highlight or treatment if available, fixed colour
    # otherwise
    {
      if (has_highlight) {
        ggplot2::geom_point(
          ggplot2::aes(colour = highlight_value),
          alpha = point_alpha,
          size  = point_size
        )
      } else if (has_treatment) {
        ggplot2::geom_point(
          ggplot2::aes(colour = treatment_group),
          alpha = point_alpha,
          size  = point_size
        )
      } else {
        ggplot2::geom_point(
          colour = point_colour,
          alpha  = point_alpha,
          size   = point_size
        )
      }
    } +

    # Colour scale: highlight (categorical or continuous) takes precedence
    # over treatment (colorblind-friendly Set2/viridis throughout)
    {
      if (has_highlight && highlight_type == "categorical") {
        n_levels <- nlevels(plot_df$highlight_value)
        if (n_levels <= 8) {
          ggplot2::scale_colour_brewer(palette = "Set2", name = colour_name)
        } else {
          ggplot2::scale_colour_viridis_d(name = colour_name)
        }
      } else if (has_highlight) {
        ggplot2::scale_colour_viridis_c(name = colour_name)
      } else if (has_treatment) {
        ggplot2::scale_colour_brewer(palette = "Set2", name = "Treatment")
      }
    } +

    # Metric annotation
    ggplot2::annotate("text",
                      x      = ann_x,
                      y      = ann_y,
                      label  = ann_label,
                      hjust  = 0,
                      vjust  = ann_vjust,
                      size   = 3.5,
                      family = "mono") +

    # Equal axis limits and scale
    ggplot2::coord_fixed(ratio = 1, xlim = lims, ylim = lims) +

    # Labels
    ggplot2::labs(
      x        = "Observed",
      y        = "CV Predicted",
      title    = "Cross-validated prediction",
      subtitle = paste0("n = ", x$n_samples, "  |  method = ",
                        x$model_params$method)
    ) +

    # Theme
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 10),
      axis.title       = ggplot2::element_text(size = 11),
      axis.text        = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = if (use_colour) "right" else "none",
      ...
    )

  p
}
