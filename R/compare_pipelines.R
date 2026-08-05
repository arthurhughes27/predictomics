# =============================================================================
# compare_pipelines.R
# Compare predictive performance across a set of pipeline configurations,
# against a baseline (covariates/mean-only) and a user-specified reference
# pipeline.
# =============================================================================


# -----------------------------------------------------------------------------
#' Compare predictive pipelines against a baseline and reference model
#'
#' @description
#' Fits a baseline model (\code{X = NULL}), a user-specified reference
#' pipeline, and a set of alternative pipelines that vary exactly one
#' pipeline stage (feature selection, feature engineering, or model choice)
#' relative to the reference, then compares their cross-validated
#' performance. Each pipeline is fit via a single sequential call to
#' \code{\link{predict_cv}}; the \code{K} alternative pipelines are always
#' fit one after another (never in parallel with each other) to avoid
#' overloading system resources, regardless of any \code{future} plan set for
#' \code{predict_cv}'s own inner CV loop.
#'
#' @details
#' Exactly one axis varies across the \code{K} alternative pipelines, selected
#' by \code{option_type}. Four axes are supported: the three pipeline stages
#' (\code{engineering_params}, \code{selection_params}, \code{model_params})
#' and the input data itself (\code{X} or \code{Y}):
#' \itemize{
#'   \item \code{"selection"}: each element of \code{option_choices}
#'     replaces \code{reference_params$selection_params}; engineering and
#'     model choice are held fixed at their reference values.
#'   \item \code{"engineering"}: each element of \code{option_choices}
#'     replaces \code{reference_params$engineering_params}; selection and
#'     model choice are held fixed.
#'   \item \code{"model"}: each element of \code{option_choices} replaces
#'     \code{reference_params$model_params}; selection and engineering are
#'     held fixed.
#'   \item \code{"predictors"}: each element of \code{option_choices} is a
#'     numeric matrix (with \code{nrow(.) == length(Y)}) that replaces
#'     \code{X} for that pipeline; \code{engineering_params},
#'     \code{selection_params}, and \code{model_params} are held fixed at
#'     their reference values, and \code{Y} is held fixed at its reference
#'     value (\code{reference_params$Y}, or the \code{Y} argument if not
#'     supplied). Useful for comparing candidate predictor panels (e.g.
#'     different omics layers, or gene subsets) under an otherwise identical
#'     pipeline.
#'   \item \code{"response"}: each element of \code{option_choices} is a
#'     numeric vector (with length matching the reference \code{X}'s
#'     \code{nrow}, or \code{length(Y)} if \code{X = NULL}) that replaces
#'     \code{Y} for that pipeline; \code{X}, \code{engineering_params},
#'     \code{selection_params}, and \code{model_params} are held fixed at
#'     their reference values (\code{reference_params$X}, or the \code{X}
#'     argument if not supplied). Useful for comparing how well the same
#'     predictor set predicts several candidate outcome variables.
#' }
#' \code{option_choices} should not include a configuration identical to the
#' reference's corresponding stage; \code{compare_pipelines()} does not check
#' for or de-duplicate this.
#'
#' The **baseline model** always uses \code{X = NULL} (see
#' \code{\link{predict_cv}}), i.e. \code{engineering_params} and
#' \code{selection_params} are both \code{NULL}, but reuses
#' \code{reference_params$model_params} as its model choice. For
#' \code{option_type = "response"}, the baseline is computed against the
#' single reference \code{Y} only (\code{reference_params$Y}, or the
#' \code{Y} argument); its metrics are therefore directly comparable to the
#' \code{"Reference"} row but not necessarily to the \code{K} response
#' options, which predict a different variable.
#'
#' **gene_level_fc row parity**: if \code{option_type = "engineering"} and
#' at least one of the reference or the \code{K} options has
#' \code{engineering_params$gene_level_fc = TRUE}, then \code{individual_id}
#' and \code{timepoint} must be supplied, and every pipeline that does
#' \strong{not} use \code{gene_level_fc} (including the baseline) has its
#' data restricted to \code{timepoint == 1} rows before fitting. This
#' matches the one-row-per-individual structure produced by
#' \code{gene_level_fc}, so that all pipelines are compared on the same
#' number of observations and CV folds are of comparable size. This
#' restriction is not applied when comparing \code{"selection"},
#' \code{"model"}, \code{"predictors"}, or \code{"response"} options
#' (\code{gene_level_fc} cannot appear as a choice there) and does not attempt
#' to reconcile with \code{selection_params$dearseq_mode = "paired"}, which
#' operates upstream of engineering and is unaffected by this restriction.
#'
#' **Reference \code{X}/\code{Y}**: for \code{option_type = "predictors"} and
#' \code{"response"} respectively, \code{reference_params$X} and
#' \code{reference_params$Y} specify the fixed value used for the
#' \code{"Reference"} row (and, for \code{"response"}, the \code{"Baseline"}
#' row). If not supplied, they default to the \code{X}/\code{Y} arguments,
#' so a call that does not set them behaves exactly as if the top-level
#' \code{X}/\code{Y} were the reference predictor set/response. These
#' elements of \code{reference_params} are ignored for all other
#' \code{option_type} values.
#'
#' **Error handling**: each of the \code{K + 2} pipeline fits (baseline,
#' reference, and the \code{K} options) is wrapped in its own error handler.
#' If a fit fails, a message is printed (regardless of \code{verbose}) and
#' that pipeline is excluded from the returned results; the remaining
#' pipelines are still attempted.
#'
#' **Message suppression**: every \code{\link{predict_cv}} call is wrapped in
#' \code{suppressMessages()}, so none of \code{predict_cv}'s own progress or
#' diagnostic messages (e.g. the double-selection note when both an explicit
#' selection method and an embedded selection model are specified) are ever
#' shown, regardless of its internal \code{verbose} gating. Only
#' \code{compare_pipelines()}'s own messages (\code{"Fitting pipeline: ..."},
#' gated by this function's \code{verbose} argument, and pipeline-failure
#' notices, always shown) are printed.
#'
#' **Same-class option labelling**: when \code{option_choices} is unnamed,
#' or contains blank/duplicate names, labels are generated from each choice's
#' \code{method} element. If multiple choices share the same method (e.g.
#' two \code{"spearman"} selections at different thresholds), they are
#' distinguished by appending \code{"_1"}, \code{"_2"}, etc, in the order
#' supplied. User-supplied names are always preferred where present, unique,
#' and non-blank.
#'
#' @param Y Numeric vector of length n. The response variable to be predicted.
#' @param X Numeric matrix of dimensions n x p. The predictor matrix used for
#'   the reference pipeline and all \code{K} options. Passed to
#'   \code{\link{predict_cv}} unchanged (the baseline model always uses
#'   \code{X = NULL} internally regardless of this argument).
#' @param option_type Character string. Which axis varies across
#'   \code{option_choices}. One of \code{"selection"}, \code{"engineering"},
#'   \code{"model"}, \code{"predictors"}, or \code{"response"}.
#' @param option_choices A list of \code{K} elements, one per alternative
#'   pipeline. Each element has the same structure as
#'   \code{selection_params} (\code{option_type = "selection"}),
#'   \code{engineering_params} (\code{"engineering"}), \code{model_params}
#'   (\code{"model"}) in \code{\link{predict_cv}}, a numeric matrix replacing
#'   \code{X} (\code{"predictors"}), or a numeric vector replacing \code{Y}
#'   (\code{"response"}). May be named to control the labels used in the
#'   returned results and plot; see Details.
#' @param reference_params A named list with elements \code{engineering_params},
#'   \code{selection_params}, \code{model_params}, \code{X}, and \code{Y},
#'   specifying the fixed reference pipeline. Whichever of these corresponds
#'   to \code{option_type} is overridden by each element of
#'   \code{option_choices} in turn; the rest are held fixed for the reference
#'   and all \code{K} options. \code{X} and \code{Y} are only relevant for
#'   \code{option_type = "predictors"}/\code{"response"} respectively (see
#'   Details) and default to the \code{X}/\code{Y} arguments when not
#'   supplied. Defaults to \code{list(engineering_params = NULL,
#'   selection_params = NULL, model_params = list(method = "lm"))}.
#' @param cv_type,folds,seed,outside_cv,treatment,treatment_predictor,covariates,individual_id,timepoint
#'   As in \code{\link{predict_cv}}, applied identically to the baseline,
#'   reference, and all \code{K} options (subject to the gene_level_fc row
#'   restriction described in Details).
#' @param metric Character string. The metric used to sort and highlight
#'   results by default in \code{\link{print.predictomics_comparison}} and
#'   \code{\link{plot.predictomics_comparison}}. One of \code{"RMSE"},
#'   \code{"sRMSE"} (default), \code{"R2"}, or \code{"SpearmanR"}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages as each
#'   pipeline is fit. Defaults to \code{TRUE}.
#'
#' @return An object of class \code{"predictomics_comparison"}, a named list
#'   containing:
#'   \describe{
#'     \item{\code{results}}{A data frame with one row per successfully-fit
#'       pipeline (\code{"Baseline"}, \code{"Reference"}, and the \code{K}
#'       option labels), and columns \code{pipeline}, \code{role} (one of
#'       \code{"baseline"}, \code{"reference"}, \code{"option"}),
#'       \code{RMSE}, \code{sRMSE}, \code{R2}, \code{SpearmanR}.}
#'     \item{\code{fits}}{A named list of the underlying \code{predictomics}
#'       objects for each successfully-fit pipeline.}
#'     \item{\code{option_type}}{The \code{option_type} argument.}
#'     \item{\code{option_choices}}{The (labelled) \code{option_choices}
#'       argument.}
#'     \item{\code{reference_params}}{The \code{reference_params} argument.}
#'     \item{\code{metric}}{The \code{metric} argument.}
#'     \item{\code{call}}{The matched call.}
#'   }
#'
#' @seealso \code{\link{predict_cv}}, \code{\link{metrics.predictomics}},
#'   \code{\link{print.predictomics_comparison}},
#'   \code{\link{plot.predictomics_comparison}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 60; p <- 100
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' colnames(X) <- paste0("gene", seq_len(p))
#' Y <- X[, 1] * 2 + rnorm(n)
#'
#' cmp <- compare_pipelines(
#'   Y = Y, X = X,
#'   option_type    = "selection",
#'   option_choices = list(
#'     list(method = "pearson",  top_n = 20),
#'     list(method = "spearman", top_n = 20),
#'     list(method = "variance", top_n = 20)
#'   ),
#'   reference_params = list(
#'     selection_params = list(method = "pearson", top_n = 50),
#'     model_params     = list(method = "glmnet")
#'   )
#' )
#' print(cmp)
#' plot(cmp)
#'
#' # Compare candidate predictor panels (X varies, Y fixed)
#' cmp_x <- compare_pipelines(
#'   Y = Y, X = X,
#'   option_type    = "predictors",
#'   option_choices = list(
#'     first_half  = X[, 1:50],
#'     second_half = X[, 51:100]
#'   ),
#'   reference_params = list(model_params = list(method = "glmnet"))
#' )
#'
#' # Compare candidate response variables (Y varies, X fixed)
#' Y2 <- X[, 2] * 2 + rnorm(n)
#' cmp_y <- compare_pipelines(
#'   Y = Y, X = X,
#'   option_type    = "response",
#'   option_choices = list(alt_response = Y2),
#'   reference_params = list(model_params = list(method = "glmnet"))
#' )
#' }
#'
#' @export
# -----------------------------------------------------------------------------
compare_pipelines <- function(Y,
                              X,
                              option_type,
                              option_choices,
                              reference_params    = list(
                                engineering_params = NULL,
                                selection_params   = NULL,
                                model_params        = list(method = "lm")
                              ),
                              cv_type             = "kfold",
                              folds               = 10L,
                              seed                = 12345L,
                              outside_cv          = FALSE,
                              treatment           = NULL,
                              treatment_predictor = FALSE,
                              covariates          = NULL,
                              individual_id       = NULL,
                              timepoint           = NULL,
                              metric              = "sRMSE",
                              verbose             = TRUE) {

  cl <- match.call()

  # ---------------------------------------------------------------------------
  # 1. Validate
  # ---------------------------------------------------------------------------
  if (!option_type %in% c("selection", "engineering", "model", "predictors",
                          "response"))
    stop("[predictomics] option_type must be one of 'selection', ",
         "'engineering', 'model', 'predictors', or 'response'.", call. = FALSE)

  if (!is.list(option_choices) || length(option_choices) == 0L)
    stop("[predictomics] option_choices must be a non-empty list.",
         call. = FALSE)

  if (!is.list(reference_params))
    stop("[predictomics] reference_params must be a named list.",
         call. = FALSE)

  if (!metric %in% c("RMSE", "sRMSE", "R2", "SpearmanR"))
    stop("[predictomics] metric must be one of 'RMSE', 'sRMSE', 'R2', or ",
         "'SpearmanR'.", call. = FALSE)

  ref_engineering <- reference_params$engineering_params
  ref_selection   <- reference_params$selection_params
  ref_model       <- reference_params$model_params %||% list(method = "lm")
  ref_X           <- reference_params$X %||% X
  ref_Y           <- reference_params$Y %||% Y

  # reference_params$X / $Y, if supplied, must be dimensionally consistent
  # with the reference Y / X respectively (defaulting to the top-level
  # arguments when not supplied, so a call that doesn't set them behaves as
  # if X/Y were themselves the reference predictor set/response)
  if (!is.null(reference_params$X)) {
    if (!is.matrix(reference_params$X) || !is.numeric(reference_params$X) ||
        nrow(reference_params$X) != length(Y))
      stop("[predictomics] reference_params$X, if supplied, must be a ",
           "numeric matrix with nrow(.) == length(Y) (", length(Y), ").",
           call. = FALSE)
  }
  if (!is.null(reference_params$Y)) {
    n_ref <- if (!is.null(X)) nrow(X) else length(Y)
    if (!is.numeric(reference_params$Y) || !is.null(dim(reference_params$Y)) ||
        length(reference_params$Y) != n_ref)
      stop("[predictomics] reference_params$Y, if supplied, must be a ",
           "numeric vector with length ", n_ref, " (nrow(X), or length(Y) ",
           "if X is NULL).", call. = FALSE)
  }

  # option_choices content validation for the data-varying option types
  if (option_type == "predictors") {
    bad_type <- !vapply(option_choices, function(o)
      is.matrix(o) && is.numeric(o), logical(1))
    if (any(bad_type))
      stop("[predictomics] For option_type = 'predictors', every element of ",
           "option_choices must be a numeric matrix.", call. = FALSE)
    bad_n <- vapply(option_choices, function(o) nrow(o) != length(Y), logical(1))
    if (any(bad_n))
      stop("[predictomics] For option_type = 'predictors', every element of ",
           "option_choices must have nrow(.) == length(Y) (", length(Y),
           ").", call. = FALSE)
  }
  if (option_type == "response") {
    bad_type <- !vapply(option_choices, function(o)
      is.numeric(o) && is.null(dim(o)), logical(1))
    if (any(bad_type))
      stop("[predictomics] For option_type = 'response', every element of ",
           "option_choices must be a numeric vector.", call. = FALSE)
    n_ref <- if (!is.null(X)) nrow(X) else length(Y)
    bad_n <- vapply(option_choices, function(o) length(o) != n_ref, logical(1))
    if (any(bad_n))
      stop("[predictomics] For option_type = 'response', every element of ",
           "option_choices must have length ", n_ref, " (nrow(X), or ",
           "length(Y) if X is NULL).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 2. Label the K options
  # ---------------------------------------------------------------------------
  fallback_label <- if (option_type %in% c("predictors", "response"))
    option_type else "option"
  option_labels  <- .make_option_labels(option_choices, fallback_label)
  names(option_choices) <- option_labels

  # ---------------------------------------------------------------------------
  # 3. gene_level_fc row-parity setup
  # ---------------------------------------------------------------------------
  option_uses_gene_level_fc <- if (option_type == "engineering") {
    vapply(option_choices, .uses_gene_level_fc, logical(1))
  } else {
    rep(FALSE, length(option_choices))
  }
  ref_uses_gene_level_fc <- option_type == "engineering" &&
    .uses_gene_level_fc(ref_engineering)

  needs_row_parity <- any(option_uses_gene_level_fc) || ref_uses_gene_level_fc

  if (needs_row_parity)
    .validate_individual_timepoint_pairing(
      individual_id, timepoint, length(Y),
      context = "compare_pipelines() with a gene_level_fc option"
    )

  # ---------------------------------------------------------------------------
  # 4. Assemble the K + 2 pipeline specifications
  # ---------------------------------------------------------------------------
  specs <- list()

  specs[["Baseline"]] <- list(
    role               = "baseline",
    engineering_params = NULL,
    selection_params   = NULL,
    model_params       = ref_model,
    X_override         = NULL,
    Y_override         = ref_Y,
    use_X              = FALSE,
    uses_gene_level_fc = FALSE
  )

  specs[["Reference"]] <- list(
    role               = "reference",
    engineering_params = ref_engineering,
    selection_params   = ref_selection,
    model_params       = ref_model,
    X_override         = ref_X,
    Y_override         = ref_Y,
    use_X              = TRUE,
    uses_gene_level_fc = ref_uses_gene_level_fc
  )

  for (i in seq_along(option_choices)) {
    lbl <- option_labels[i]
    specs[[lbl]] <- switch(
      option_type,
      selection = list(
        role               = "option",
        engineering_params = ref_engineering,
        selection_params   = option_choices[[i]],
        model_params       = ref_model,
        X_override         = NULL,
        Y_override         = NULL,
        use_X              = TRUE,
        uses_gene_level_fc = FALSE
      ),
      engineering = list(
        role               = "option",
        engineering_params = option_choices[[i]],
        selection_params   = ref_selection,
        model_params       = ref_model,
        X_override         = NULL,
        Y_override         = NULL,
        use_X              = TRUE,
        uses_gene_level_fc = option_uses_gene_level_fc[i]
      ),
      model = list(
        role               = "option",
        engineering_params = ref_engineering,
        selection_params   = ref_selection,
        model_params       = option_choices[[i]],
        X_override         = NULL,
        Y_override         = NULL,
        use_X              = TRUE,
        uses_gene_level_fc = FALSE
      ),
      predictors = list(
        role               = "option",
        engineering_params = ref_engineering,
        selection_params   = ref_selection,
        model_params       = ref_model,
        X_override         = option_choices[[i]],
        Y_override         = NULL,
        use_X              = TRUE,
        uses_gene_level_fc = FALSE
      ),
      response = list(
        role               = "option",
        engineering_params = ref_engineering,
        selection_params   = ref_selection,
        model_params       = ref_model,
        X_override         = NULL,
        Y_override         = option_choices[[i]],
        use_X              = TRUE,
        uses_gene_level_fc = FALSE
      )
    )
  }

  # ---------------------------------------------------------------------------
  # 5. Fit each pipeline sequentially, catching and continuing on error
  # ---------------------------------------------------------------------------
  fits    <- list()
  rows    <- list()

  for (lbl in names(specs)) {

    spec <- specs[[lbl]]

    if (verbose)
      message("[predictomics] Fitting pipeline: ", lbl, " (", spec$role, ")")

    fit <- tryCatch({

      restrict <- needs_row_parity && !spec$uses_gene_level_fc

      this_X <- spec$X_override %||% X
      this_Y <- spec$Y_override %||% Y

      data <- if (restrict) {
        .restrict_to_timepoint1(this_Y, this_X, covariates, treatment, timepoint)
      } else {
        list(Y = this_Y, X = this_X, covariates = covariates, treatment = treatment)
      }

      pass_individual_id <- if (restrict) NULL else individual_id
      pass_timepoint      <- if (restrict) NULL else timepoint

      suppressMessages(predict_cv(
        Y                   = data$Y,
        X                   = if (spec$use_X) data$X else NULL,
        cv_type             = cv_type,
        folds               = folds,
        seed                = seed,
        engineering_params  = spec$engineering_params,
        selection_params    = spec$selection_params,
        model_params        = spec$model_params,
        outside_cv          = outside_cv,
        treatment           = data$treatment,
        treatment_predictor = treatment_predictor,
        covariates          = data$covariates,
        individual_id       = pass_individual_id,
        timepoint           = pass_timepoint,
        verbose             = FALSE
      ))

    }, error = function(e) {
      message("[predictomics] Pipeline '", lbl, "' failed and was excluded ",
              "from the results: ", conditionMessage(e))
      NULL
    })

    if (!is.null(fit)) {
      fits[[lbl]] <- fit
      m <- metrics.predictomics(fit, digits = 4)
      rows[[lbl]] <- data.frame(
        pipeline  = lbl,
        role      = spec$role,
        RMSE      = unname(m["RMSE"]),
        sRMSE     = unname(m["sRMSE"]),
        R2        = unname(m["R2"]),
        SpearmanR = unname(m["SpearmanR"]),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L)
    stop("[predictomics] All pipelines failed; no results to report.",
         call. = FALSE)

  results <- do.call(rbind, rows)
  rownames(results) <- NULL

  # ---------------------------------------------------------------------------
  # 6. Return
  # ---------------------------------------------------------------------------
  structure(
    list(
      results          = results,
      fits             = fits,
      option_type      = option_type,
      option_choices   = option_choices,
      reference_params = reference_params,
      metric           = metric,
      call             = cl
    ),
    class = "predictomics_comparison"
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# -----------------------------------------------------------------------------
#' Detect whether an engineering_params list uses gene_level_fc
#'
#' @param engineering_params An \code{engineering_params} list, or \code{NULL}.
#' @return Logical scalar.
#' @keywords internal
# -----------------------------------------------------------------------------
.uses_gene_level_fc <- function(engineering_params) {
  isTRUE(engineering_params$gene_level_fc)
}


# -----------------------------------------------------------------------------
#' Human-readable label for an option_type value
#'
#' @description
#' Maps \code{option_type} to the phrase used in
#' \code{\link{plot.predictomics_comparison}}'s title and
#' \code{\link{print.predictomics_comparison}}'s summary.
#'
#' @param option_type Character string. One of \code{"selection"},
#'   \code{"engineering"}, \code{"model"}, \code{"predictors"}, or
#'   \code{"response"}.
#' @return A character string.
#' @keywords internal
# -----------------------------------------------------------------------------
.option_type_label <- function(option_type) {
  switch(
    option_type,
    selection   = "feature selection",
    engineering = "engineering",
    model       = "model choice",
    predictors  = "predictor set",
    response    = "response variable"
  )
}


# -----------------------------------------------------------------------------
#' Restrict Y, X, covariates, and treatment to post-treatment rows
#'
#' @description
#' Subsets to \code{timepoint == 1} rows, used by \code{\link{compare_pipelines}}
#' to keep non-\code{gene_level_fc} pipelines comparable (in row count and CV
#' fold structure) to a \code{gene_level_fc} pipeline, which collapses paired
#' rows to one row per individual.
#'
#' @param Y Numeric response vector.
#' @param X Numeric predictor matrix, or \code{NULL}.
#' @param covariates A covariate matrix or data frame, or \code{NULL}.
#' @param treatment A treatment vector, or \code{NULL}.
#' @param timepoint A binary numeric vector (0/1), paired with \code{Y}.
#' @return A named list with the restricted \code{Y}, \code{X},
#'   \code{covariates}, and \code{treatment}.
#' @keywords internal
# -----------------------------------------------------------------------------
.restrict_to_timepoint1 <- function(Y, X, covariates, treatment, timepoint) {

  keep <- timepoint == 1

  list(
    Y          = Y[keep],
    X          = if (!is.null(X)) X[keep, , drop = FALSE] else NULL,
    covariates = if (!is.null(covariates)) {
      if (is.data.frame(covariates)) covariates[keep, , drop = FALSE]
      else covariates[keep, , drop = FALSE]
    } else NULL,
    treatment  = if (!is.null(treatment)) treatment[keep] else NULL
  )
}


# -----------------------------------------------------------------------------
#' Generate distinguishing labels for a list of option_choices
#'
#' @description
#' Resolves display labels for each element of \code{option_choices} passed
#' to \code{\link{compare_pipelines}}. User-supplied names are used where
#' present, unique, and non-blank. Remaining elements fall back to their
#' \code{method} element (for parameter-list options) or to
#' \code{fallback_label} (for matrix/vector options, i.e.
#' \code{option_type \%in\% c("predictors", "response")}); if this produces
#' duplicates (e.g. two choices with \code{method = "spearman"}, or two
#' unnamed predictor-set options), duplicates are disambiguated by appending
#' \code{"_1"}, \code{"_2"}, etc, in order of appearance. As a final
#' safeguard, any labels still duplicated after this (e.g. a user-supplied
#' name colliding with a generated one) have their index appended.
#'
#' @param option_choices A list of parameter lists, matrices, or vectors,
#'   optionally named.
#' @param fallback_label Character string used as the base label for
#'   elements without a \code{method} element (or that are not lists at
#'   all). Defaults to \code{"option"}.
#' @return A character vector of labels, length \code{length(option_choices)}.
#' @keywords internal
# -----------------------------------------------------------------------------
.make_option_labels <- function(option_choices, fallback_label = "option") {

  n <- length(option_choices)
  user_names <- names(option_choices)
  if (is.null(user_names)) user_names <- rep("", n)
  user_names[is.na(user_names)] <- ""

  methods <- vapply(option_choices, function(p) {
    m <- if (is.list(p)) p$method else NULL
    if (is.null(m) || !is.character(m) || length(m) != 1L || !nzchar(m))
      fallback_label
    else
      m
  }, character(1))

  has_user_name <- nzchar(user_names) & !duplicated(user_names) &
    !duplicated(user_names, fromLast = TRUE)
  # A blank name never counts as a valid user name
  has_user_name[!nzchar(user_names)] <- FALSE

  labels <- ifelse(has_user_name, user_names, methods)

  # Method-grouped numbering for the auto-generated (non-user-named) labels
  dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(dup & !has_user_name)) {
    counters <- list()
    for (i in seq_len(n)) {
      if (dup[i] && !has_user_name[i]) {
        key <- methods[i]
        counters[[key]] <- (counters[[key]] %||% 0L) + 1L
        labels[i] <- paste0(key, "_", counters[[key]])
      }
    }
  }

  # Final tiebreaker: append index to any labels still duplicated
  still_dup <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  if (any(still_dup))
    labels[still_dup] <- paste0(labels[still_dup], "_", seq_len(n)[still_dup])

  labels
}


# -----------------------------------------------------------------------------
#' Print a predictomics_comparison object
#'
#' @description
#' Prints a concise summary table comparing the baseline, reference, and
#' option pipelines fit by \code{\link{compare_pipelines}}.
#'
#' @param x A \code{predictomics_comparison} object returned by
#'   \code{\link{compare_pipelines}}.
#' @param digits Integer. Number of decimal places for metrics. Defaults to
#'   \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{compare_pipelines}},
#'   \code{\link{plot.predictomics_comparison}}
#'
#' @export
# -----------------------------------------------------------------------------
print.predictomics_comparison <- function(x, digits = 4, ...) {

  if (!inherits(x, "predictomics_comparison"))
    stop("[predictomics] x must be a predictomics_comparison object returned ",
         "by compare_pipelines().", call. = FALSE)

  res <- x$results
  res[c("RMSE", "sRMSE", "R2", "SpearmanR")] <- lapply(
    res[c("RMSE", "sRMSE", "R2", "SpearmanR")], round, digits = digits
  )
  res <- res[order(res[[x$metric]]), ]

  cat("\n")
  cat("=================================================\n")
  cat(" predictomics: pipeline comparison result \n")
  cat("=================================================\n")

  cat("\nCall:\n ")
  cat(deparse(x$call), "\n")

  cat("\nOption type :", x$option_type, "\n")
  cat("Options     :", paste(names(x$option_choices), collapse = ", "), "\n")
  cat("Sorted by   :", x$metric, "\n\n")

  print(res[, c("pipeline", "role", "RMSE", "sRMSE", "R2", "SpearmanR")],
        row.names = FALSE)

  cat("\n=================================================\n\n")

  invisible(x)
}


# -----------------------------------------------------------------------------
#' Plot a predictomics_comparison object
#'
#' @description
#' Produces a grouped bar chart comparing the reference and alternative
#' pipelines fit by \code{\link{compare_pipelines}} on one or all performance
#' metrics, against the baseline shown as a dashed reference line.
#'
#' @details
#' The Reference and Alternative pipelines are shown as coloured bars; the
#' Baseline is shown instead as a dashed horizontal line (both are part of a
#' single "Pipelines" legend) annotated with its metric value at the right of
#' the panel. The best-performing bar (lowest \code{RMSE}/\code{sRMSE}, or
#' highest \code{R2}/\code{SpearmanR}) is marked with a diagonal hatch pattern
#' (via \pkg{ggpattern}), with its own "Best pipeline" legend key; the
#' "Pipelines" (fill) legend keys never show hatching themselves, regardless
#' of which bar is best, so the two legends stay visually distinct. When
#' \code{metric = "all"}, the four metrics are shown as facets, each with its
#' own best-bar highlight and baseline annotation.
#'
#' @param x A \code{predictomics_comparison} object returned by
#'   \code{\link{compare_pipelines}}.
#' @param metric Character string. Which metric to plot. One of
#'   \code{"RMSE"}, \code{"sRMSE"}, \code{"R2"}, \code{"SpearmanR"}, or
#'   \code{"all"} (facet over all four). Defaults to \code{x$metric}.
#' @param sort Logical. If \code{TRUE}, pipelines are ordered by the plotted
#'   metric (or by \code{x$metric} when \code{metric = "all"}). Defaults to
#'   \code{FALSE}, keeping the order in which pipelines were supplied
#'   (Reference, then the \code{K} options in order; the Baseline is always
#'   shown as a line rather than a bar).
#' @param ... Additional arguments passed to \code{ggplot2::theme}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{compare_pipelines}},
#'   \code{\link{print.predictomics_comparison}}
#'
#' @export
# -----------------------------------------------------------------------------
plot.predictomics_comparison <- function(x,
                                         metric = x$metric,
                                         sort   = FALSE,
                                         ...) {

  if (!inherits(x, "predictomics_comparison"))
    stop("[predictomics] x must be a predictomics_comparison object returned ",
         "by compare_pipelines().", call. = FALSE)

  if (!metric %in% c("RMSE", "sRMSE", "R2", "SpearmanR", "all"))
    stop("[predictomics] metric must be one of 'RMSE', 'sRMSE', 'R2', ",
         "'SpearmanR', or 'all'.", call. = FALSE)

  res <- x$results
  res$role <- factor(res$role, levels = c("baseline", "reference", "option"),
                     labels = c("Baseline", "Reference", "Alternative"))

  sort_metric <- if (metric == "all") x$metric else metric

  if (sort) {
    ord <- order(res[[sort_metric]])
    res$pipeline <- factor(res$pipeline, levels = res$pipeline[ord])
  } else {
    res$pipeline <- factor(res$pipeline, levels = res$pipeline)
  }

  baseline_row <- res[res$role == "Baseline", , drop = FALSE]
  bars_res     <- res[res$role != "Baseline", , drop = FALSE]
  bars_res$role <- droplevels(bars_res$role)
  bars_res$pipeline <- droplevels(bars_res$pipeline)

  option_label <- .option_type_label(x$option_type)
  title_text   <- paste0("Pipeline comparison: ", option_label)

  pipeline_colours <- c(Reference = "#FC8D62", Alternative = "#8DA0CB")
  bar_width        <- 0.7

  # Bar geom: a hatched pattern marks the best-performing bar. The "pattern"
  # legend shows only the "Best pipeline" key; the "Pipelines" (fill) legend
  # keys are overridden to never show hatching, regardless of is_best, so the
  # two legends don't visually interfere with each other.
  bar_geom <- ggpattern::geom_col_pattern(
    ggplot2::aes(pattern = is_best),
    alpha            = 0.9, width = bar_width,
    pattern_fill     = "grey20", pattern_colour = NA,
    pattern_density  = 0.25, pattern_spacing = 0.03, pattern_angle = 45,
    pattern_key_scale_factor = 0.6
  )

  pattern_scale <- ggpattern::scale_pattern_manual(
    name   = NULL,
    values = c(`TRUE` = "stripe", `FALSE` = "none"),
    breaks = "TRUE",
    labels = "Best pipeline",
    guide  = ggplot2::guide_legend(
      override.aes = list(fill = "grey90", pattern_fill = "grey20"),
      order        = 2
    )
  )

  fill_guide <- ggplot2::guide_legend(
    override.aes = list(pattern = "none"), order = 1
  )

  if (metric == "all") {

    plot_df <- tidyr::pivot_longer(
      bars_res,
      cols      = c("RMSE", "sRMSE", "R2", "SpearmanR"),
      names_to  = "metric",
      values_to = "value"
    )
    plot_df$metric <- factor(plot_df$metric,
                             levels = c("RMSE", "sRMSE", "R2", "SpearmanR"))

    ref_lines <- data.frame(
      metric = factor(c("RMSE", "sRMSE", "R2", "SpearmanR"),
                      levels = c("RMSE", "sRMSE", "R2", "SpearmanR")),
      value  = c(baseline_row$RMSE, baseline_row$sRMSE,
                baseline_row$R2, baseline_row$SpearmanR)
    )
    ref_lines$label <- paste0("Baseline: ", round(ref_lines$value, 3))

    best_df <- do.call(rbind, lapply(levels(plot_df$metric), function(m) {
      sub           <- plot_df[plot_df$metric == m, , drop = FALSE]
      higher_better <- m %in% c("R2", "SpearmanR")
      sub[which.max(if (higher_better) sub$value else -sub$value), , drop = FALSE]
    }))

    plot_df$is_best <- paste(plot_df$metric, plot_df$pipeline) %in%
      paste(best_df$metric, best_df$pipeline)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pipeline, y = value,
                                               fill = role)) +
      bar_geom +
      ggplot2::geom_hline(
        data        = ref_lines,
        ggplot2::aes(yintercept = value, linetype = "Baseline"),
        colour      = "grey40", linewidth = 0.6, inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data        = ref_lines,
        ggplot2::aes(x = Inf, y = value, label = label),
        hjust = 1.05, vjust = -0.5, size = 2.8, colour = "grey40",
        inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(~metric, scales = "free_y") +
      ggplot2::scale_linetype_manual(name = NULL,
                                     values = c(Baseline = "dashed")) +
      pattern_scale +
      ggplot2::labs(
        x        = "Pipeline specification", y = NULL, fill = "Pipelines",
        title    = title_text
      )

  } else {

    plot_df <- bars_res
    plot_df$metric_value <- bars_res[[metric]]

    higher_better <- metric %in% c("R2", "SpearmanR")
    best_df <- plot_df[which.max(if (higher_better) plot_df$metric_value
                                 else -plot_df$metric_value), , drop = FALSE]
    plot_df$is_best <- plot_df$pipeline %in% best_df$pipeline

    baseline_value <- baseline_row[[metric]]
    baseline_label <- paste0("Baseline: ", round(baseline_value, 3))
    baseline_df    <- data.frame(value = baseline_value)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pipeline, y = metric_value,
                                               fill = role)) +
      bar_geom +
      ggplot2::geom_hline(
        data        = baseline_df,
        ggplot2::aes(yintercept = value, linetype = "Baseline"),
        colour      = "grey40", linewidth = 0.6, inherit.aes = FALSE
      ) +
      ggplot2::annotate(
        "text", x = Inf, y = baseline_value, label = baseline_label,
        hjust = 1.05, vjust = -0.6, size = 3, colour = "grey40"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = round(metric_value, 3)),
        vjust = -0.4, size = 3, colour = "grey20"
      ) +
      ggplot2::scale_linetype_manual(name = NULL,
                                     values = c(Baseline = "dashed")) +
      pattern_scale +
      ggplot2::labs(
        x        = "Pipeline specification", y = metric, fill = "Pipelines",
        title    = title_text,
        subtitle = paste0("metric: ", metric)
      )
  }

  p <- p +
    ggplot2::scale_fill_manual(values = pipeline_colours, guide = fill_guide) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 12),
      axis.title       = ggplot2::element_text(size = 13, face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y      = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      ...
    )

  p
}
