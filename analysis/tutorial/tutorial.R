# =============================================================================
# tutorial.R
#
# A self-contained tutorial demonstrating the two main entry points of the
# predictomics package - predict_cv() and compare_pipelines() - on simulated
# data, mirroring the style of the SurrogateRank tutorial (chapter 3).
#
# Run this script top to bottom from the package root (or set your working
# directory to analysis/tutorial/ and adjust FIGURES_DIR/OUTPUTS_DIR below).
# It writes:
#   - figures/*.pdf   : one PDF per plot, ready for \includegraphics{}
#   - outputs/*.txt   : one plain-text file per console output, ready to
#                        paste into a \begin{minted}{text} ... \end{minted}
#                        block (or \lstinputlisting{}/\input{} if preferred)
#
# See README.md in this folder for guidance on lifting code/output/figures
# from this script directly into the thesis chapter.
# =============================================================================

library(predictomics)

FIGURES_DIR <- "analysis/tutorial/figures"
OUTPUTS_DIR <- "analysis/tutorial/outputs"
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Small helpers: save a plot to figures/, save captured console text to
# outputs/. Keeping these separate from the "narrative" code below keeps the
# demonstration code itself close to what you'd actually paste into a
# \begin{minted}{r} block in the chapter.
save_figure <- function(plot_obj, file, width = 7, height = 5) {
  ggplot2::ggsave(
    filename = file.path(FIGURES_DIR, file),
    plot     = plot_obj,
    width    = width,
    height   = height,
    device   = "pdf"
  )
}

save_console <- function(text, file) {
  writeLines(text, file.path(OUTPUTS_DIR, file))
  invisible(text)
}


# =============================================================================
# Part 1: predict_cv() - a single cross-validated pipeline
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Simulate data
#
# simulate_predictomics_data() generates gene expression data (X), a
# continuous response (Y), a binary treatment indicator, covariates (age,
# sex), and a partition of the genes into named genesets. Two of the ten
# genesets are "signal" genesets whose mean expression drives Y; the rest
# are pure noise, so we know in advance which features a well-behaved
# pipeline should recover.
# -----------------------------------------------------------------------------
sim <- simulate_predictomics_data(
  n                 = 60,  # samples
  p                 = 200, # genes
  n_genesets        = 10,
  geneset_size      = 20,
  n_signal_genesets = 2,
  seed              = 1
)

save_console(capture.output(str(sim)), "01_simulate_str.txt")

# -----------------------------------------------------------------------------
# 1.2 Fit a single pipeline with predict_cv()
#
# We aggregate genes into geneset-level mean expression (z-scored beforehand)
# and fit a lasso model, which performs embedded feature selection. Both the
# engineering and model steps are refit within every CV fold by default
# (outside_cv = FALSE), so the reported performance is not optimistically
# biased. model_params$scale = TRUE and compute_importance = TRUE additionally
# let us interpret which genesets the model actually relied on.
# -----------------------------------------------------------------------------
engineering_params <- list(
  method        = "engineer",
  col_transform = "z",
  genesets      = sim$genesets,
  agg_method    = "mean"
)

model_params <- list(
  method             = "lasso",
  scale              = TRUE,
  compute_importance = TRUE
)

fit <- predict_cv(
  Y                  = sim$Y,
  X                  = sim$X,
  engineering_params = engineering_params,
  model_params       = model_params,
  covariates         = sim$covariates,
  cv_type            = "kfold",
  folds              = 10,
  seed               = 12345,
  verbose            = FALSE
)

# -----------------------------------------------------------------------------
# 1.3 Inspect the result
# -----------------------------------------------------------------------------
save_console(capture.output(print(fit)), "02_predict_cv_print.txt")

p_fit <- plot(fit)
save_figure(p_fit, "tutorial_predict_cv_scatter.pdf")

# Embedded selection stability: how consistently the lasso picked out each
# geneset across the 10 outer folds.
p_stability <- plot_selection_stability(fit, type = "embedded")
save_figure(p_stability$frequency, "tutorial_selection_stability.pdf")

# Feature importance: mean standardised |coefficient| across folds, which is
# meaningful regardless of how many genesets ended up being selected.
p_importance <- plot_feature_importance(fit)
save_figure(p_importance, "tutorial_feature_importance.pdf")

# The genesets that actually generated the signal, for comparison against
# the plots above.
save_console(capture.output(print(sim$signal_genesets)), "03_signal_genesets.txt")


# =============================================================================
# Part 2: compare_pipelines() - comparing alternative pipelines
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 Compare three modelling choices against the lasso reference
#
# compare_pipelines() refits a reference pipeline plus one pipeline per entry
# of option_choices, all under an identical CV split, and returns a results
# table plus a comparison plot. Here we vary option_type = "model", holding
# the geneset-level engineering step from Part 1 fixed for every pipeline.
# -----------------------------------------------------------------------------
cmp <- compare_pipelines(
  Y                  = sim$Y,
  X                  = sim$X,
  option_type        = "model",
  option_choices     = list(
    "Random forest" = list(method = "ranger"),
    "SVR"           = list(method = "svr", scale = TRUE)
  ),
  reference_params   = list(
    engineering_params = engineering_params,
    selection_params   = NULL,
    model_params        = model_params
  ),
  covariates         = sim$covariates,
  cv_type            = "kfold",
  folds              = 10,
  seed               = 12345,
  metric             = "sRMSE",
  verbose            = FALSE
)

# -----------------------------------------------------------------------------
# 2.2 Inspect the result
# -----------------------------------------------------------------------------
save_console(capture.output(print(cmp)), "04_compare_pipelines_print.txt")

# A preview of the full results table, one row per pipeline (baseline,
# reference, and each option).
save_console(capture.output(print(cmp$results, row.names = FALSE)),
            "05_compare_pipelines_results.txt")

p_cmp <- plot(cmp)
save_figure(p_cmp, "tutorial_compare_pipelines.pdf", width = 8, height = 5)

message("Tutorial complete. Figures written to ", FIGURES_DIR,
       "; console outputs written to ", OUTPUTS_DIR, ".")
