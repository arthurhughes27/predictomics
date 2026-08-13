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
#' **Paired row-discard parity**: some pipeline configurations discard
#' pre-treatment (\code{timepoint == 0}) rows internally, either by
#' collapsing paired rows into one fold-change row per individual
#' (\code{engineering_params$gene_level_fc = TRUE}) or by using both arms for
#' an upfront filtration step and then modelling on post-treatment rows only
#' (\code{selection_params$rise_paired = TRUE}, or
#' \code{selection_params$method = "dearseq"} with
#' \code{dearseq_mode = "paired"}; see \code{\link{predict_cv}}). This applies
#' to the reference and every \code{K} option's \strong{effective}
#' \code{engineering_params}/\code{selection_params} - i.e. whichever of
#' \code{option_choices} or the corresponding \code{reference_params} element
#' is in force for that pipeline - regardless of \code{option_type}. If
#' \strong{any} pipeline (reference, baseline, or an option) uses one of
#' these, \code{individual_id} and \code{timepoint} must be supplied, and
#' every pipeline that does \strong{not} itself discard pre-treatment rows
#' (including the baseline) has its data restricted to \code{timepoint == 1}
#' rows before fitting; pipelines that do discard rows internally receive the
#' full, unrestricted data instead, so their own paired filtration step still
#' sees both arms. This keeps every pipeline being compared on the same final
#' number of observations, with comparable CV folds, even when
#' \code{option_choices} mixes paired and non-paired configurations (e.g.
#' comparing \code{selection_params$method = "dearseq"} with
#' \code{dearseq_mode = "paired"} against a plain correlation filter).
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
#' **Memory usage**: fitting many pipelines against a large \code{X} (many
#' features and/or many outer folds) can use substantial memory, since every
#' successfully-fit pipeline's full \code{predictomics} object is retained in
#' the returned \code{fits} list for the lifetime of the call. Two things
#' help: \code{diagnostics = "summary"} (see \code{diagnostics} below) drops
#' the largest per-fit fields, and \code{gc()} is called after each of the
#' \code{K + 2} pipeline fits so a finished pipeline's temporaries (its
#' engineered/restricted data copy, and, in \code{"full"} mode, its
#' untrimmed diagnostics) are reclaimed before the next, potentially equally
#' large, pipeline is fit, rather than accumulating across the sequential
#' loop. If the inner CV loop itself is parallelised (see
#' \code{\link{predict_cv}}'s Parallelisation section), note that each
#' worker process holds its own copy of \code{X} for the duration of that
#' \code{predict_cv} call; a large parallel worker count multiplies memory
#' use accordingly and is independent of the settings here.
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
#'   reference, and all \code{K} options (subject to the paired row-discard
#'   parity restriction described in Details).
#' @param metric Character string. The metric used to sort and highlight
#'   results by default in \code{\link{print.predictomics_comparison}} and
#'   \code{\link{plot.predictomics_comparison}}. One of \code{"RMSE"},
#'   \code{"sRMSE"} (default), \code{"R2"}, or \code{"SpearmanR"}.
#' @param diagnostics Character string. One of \code{"full"} (default) or
#'   \code{"summary"}. In \code{"summary"} mode, the full-length
#'   \code{selection_scores} vector (one score per candidate feature or
#'   geneset, not just the selected ones) is dropped from
#'   \code{dearseq_selection}, each element of
#'   \code{fold_selection_diagnostics}, and \code{outside_cv_selection} in
#'   every stored fit before it is added to the returned \code{fits} list;
#'   \code{selected_features} and \code{n_selected} are always kept.
#'   \code{results} (the summary metrics table) is unaffected either way,
#'   since it never depends on these fields. See Details ("Memory usage").
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
#'       objects for each successfully-fit pipeline. See \code{diagnostics}
#'       above for the \code{"summary"} mode that trims their largest
#'       fields.}
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
#'
#' # Compare selection methods where some are paired (discard pre-treatment
#' # rows internally) and some are not - individual_id/timepoint only need to
#' # be supplied once; row-parity restriction is handled automatically.
#' n_ind <- 30
#' individual_id <- rep(seq_len(n_ind), each = 2)
#' timepoint     <- rep(c(0, 1), times = n_ind)
#' cmp_paired <- compare_pipelines(
#'   Y = Y[seq_len(n_ind * 2)], X = X[seq_len(n_ind * 2), ],
#'   option_type    = "selection",
#'   option_choices = list(
#'     dearseq_paired = list(method = "dearseq", dearseq_mode = "paired",
#'                          threshold = 0.05),
#'     correlation     = list(method = "pearson", top_n = 20)
#'   ),
#'   reference_params = list(
#'     selection_params = list(method = "variance", top_n = 20),
#'     model_params      = list(method = "lm")
#'   ),
#'   individual_id = individual_id,
#'   timepoint     = timepoint
#' )
#'
#' # For a large X and/or many options, reduce memory usage by dropping the
#' # full-length per-feature selection scores from each stored fit
#' cmp_lean <- compare_pipelines(
#'   Y = Y, X = X,
#'   option_type    = "selection",
#'   option_choices = list(
#'     list(method = "pearson",  top_n = 20),
#'     list(method = "spearman", top_n = 20)
#'   ),
#'   reference_params = list(
#'     selection_params = list(method = "variance", top_n = 20),
#'     model_params      = list(method = "glmnet")
#'   ),
#'   diagnostics = "summary"
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
                              diagnostics         = "full",
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

  if (!diagnostics %in% c("full", "summary"))
    stop("[predictomics] diagnostics must be one of 'full' or 'summary'.",
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
  # 3. Paired row-discard parity setup. For each spec (reference and every
  # option), determine its *effective* engineering_params/selection_params -
  # whichever of option_choices[[i]] or the reference value is in force for
  # that spec, given option_type - and check whether that combination
  # discards pre-treatment rows internally (gene_level_fc, rise_paired, or
  # dearseq_mode = "paired"; see .discards_pretreatment_rows()). This check is
  # option_type-agnostic: e.g. for option_type = "selection", each option's
  # own selection_params is checked individually, so a mix of paired and
  # non-paired selection methods across option_choices is detected correctly.
  # ---------------------------------------------------------------------------
  option_discards_rows <- vapply(seq_along(option_choices), function(i) {
    eng <- if (option_type == "engineering") option_choices[[i]] else ref_engineering
    sel <- if (option_type == "selection")   option_choices[[i]] else ref_selection
    .discards_pretreatment_rows(eng, sel)
  }, logical(1))

  ref_discards_rows <- .discards_pretreatment_rows(ref_engineering, ref_selection)

  needs_row_parity <- any(option_discards_rows) || ref_discards_rows

  if (needs_row_parity)
    .validate_individual_timepoint_pairing(
      individual_id, timepoint, length(Y),
      context = paste0(
        "compare_pipelines() with a paired option (engineering_params$",
        "gene_level_fc, selection_params$rise_paired, or ",
        "selection_params$dearseq_mode = 'paired')"
      )
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
    discards_rows      = FALSE
  )

  specs[["Reference"]] <- list(
    role               = "reference",
    engineering_params = ref_engineering,
    selection_params   = ref_selection,
    model_params       = ref_model,
    X_override         = ref_X,
    Y_override         = ref_Y,
    use_X              = TRUE,
    discards_rows      = ref_discards_rows
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
        discards_rows      = option_discards_rows[i]
      ),
      engineering = list(
        role               = "option",
        engineering_params = option_choices[[i]],
        selection_params   = ref_selection,
        model_params       = ref_model,
        X_override         = NULL,
        Y_override         = NULL,
        use_X              = TRUE,
        discards_rows      = option_discards_rows[i]
      ),
      model = list(
        role               = "option",
        engineering_params = ref_engineering,
        selection_params   = ref_selection,
        model_params       = option_choices[[i]],
        X_override         = NULL,
        Y_override         = NULL,
        use_X              = TRUE,
        discards_rows      = option_discards_rows[i]
      ),
      predictors = list(
        role               = "option",
        engineering_params = ref_engineering,
        selection_params   = ref_selection,
        model_params       = ref_model,
        X_override         = option_choices[[i]],
        Y_override         = NULL,
        use_X              = TRUE,
        discards_rows      = option_discards_rows[i]
      ),
      response = list(
        role               = "option",
        engineering_params = ref_engineering,
        selection_params   = ref_selection,
        model_params       = ref_model,
        X_override         = NULL,
        Y_override         = option_choices[[i]],
        use_X              = TRUE,
        discards_rows      = option_discards_rows[i]
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

      restrict <- needs_row_parity && !spec$discards_rows

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
      m <- metrics.predictomics(fit, digits = 4)
      if (diagnostics == "summary") fit <- .trim_predictomics_diagnostics(fit)
      fits[[lbl]] <- fit
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

    # Each pipeline's temporaries (the restricted/engineered data copy built
    # above, and - unless diagnostics = "summary" - the untrimmed fit) can be
    # sizeable for large X; force reclamation before the next (potentially
    # equally large) pipeline is fit, rather than relying on R's lazy
    # collector to catch up across several sequential large allocations.
    gc(full = FALSE)
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
#' Detect whether a pipeline configuration discards pre-treatment rows
#'
#' @description
#' Returns \code{TRUE} if the given \code{engineering_params}/
#' \code{selection_params} combination discards pre-treatment
#' (\code{timepoint == 0}) rows internally before or during modelling - either
#' by collapsing paired rows into one fold-change row per individual
#' (\code{engineering_params$gene_level_fc = TRUE}), or by using both arms for
#' an upfront paired filtration step and then modelling on post-treatment
#' rows only (\code{selection_params$rise_paired = TRUE}, or
#' \code{selection_params$method = "dearseq"} with
#' \code{selection_params$dearseq_mode = "paired"}; see
#' \code{\link{predict_cv}}). Used by \code{\link{compare_pipelines}} to
#' detect when cross-pipeline row-parity restriction is needed.
#'
#' @param engineering_params An \code{engineering_params} list, or \code{NULL}.
#' @param selection_params A \code{selection_params} list, or \code{NULL}.
#' @return Logical scalar.
#' @keywords internal
# -----------------------------------------------------------------------------
.discards_pretreatment_rows <- function(engineering_params, selection_params) {

  uses_gene_level_fc <- isTRUE(engineering_params$gene_level_fc)

  uses_rise_paired <- !is.null(selection_params) &&
    isTRUE(selection_params$method == "rise") &&
    isTRUE(selection_params$rise_paired)

  uses_dearseq_paired <- !is.null(selection_params) &&
    isTRUE(selection_params$method == "dearseq") &&
    identical(selection_params$dearseq_mode %||% "classic", "paired")

  uses_gene_level_fc || uses_rise_paired || uses_dearseq_paired
}


# -----------------------------------------------------------------------------
#' Strip full-length selection score vectors from a predictomics fit
#'
#' @description
#' Used by \code{\link{compare_pipelines}} when \code{diagnostics = "summary"}
#' to reduce the memory footprint of each stored fit. \code{run_selection()}
#' always returns a score for \strong{every} candidate feature (or geneset),
#' not just the selected ones; \code{\link{predict_cv}} retains one such
#' vector per outer fold in \code{fold_selection_diagnostics}, plus one more
#' in \code{dearseq_selection}/\code{outside_cv_selection} where applicable.
#' For a large \code{p} and many folds/pipelines, these accumulate quickly.
#' This strips just the \code{selection_scores} element from each of those
#' three fields (setting it to \code{NULL}), leaving \code{selected_features}
#' and \code{n_selected} - the fields \code{\link{compare_pipelines}}'s own
#' results table and \code{plot()}/\code{print()} methods rely on - intact.
#' \code{fold_embedded_selection_diagnostics} (lasso/glmnet non-zero
#' coefficients) is left untouched, since it is already limited to the
#' selected features rather than every candidate.
#'
#' @param fit A \code{predictomics} object returned by
#'   \code{\link{predict_cv}}.
#' @return The same object with \code{selection_scores} removed from
#'   \code{dearseq_selection}, each element of
#'   \code{fold_selection_diagnostics}, and \code{outside_cv_selection}.
#' @keywords internal
# -----------------------------------------------------------------------------
.trim_predictomics_diagnostics <- function(fit) {

  if (!is.null(fit$dearseq_selection))
    fit$dearseq_selection$selection_scores <- NULL

  if (!is.null(fit$fold_selection_diagnostics))
    fit$fold_selection_diagnostics <- lapply(
      fit$fold_selection_diagnostics,
      function(d) {
        if (!is.null(d)) d$selection_scores <- NULL
        d
      }
    )

  if (!is.null(fit$outside_cv_selection))
    fit$outside_cv_selection$selection_scores <- NULL

  fit
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
#' to keep pipelines that do not themselves discard pre-treatment rows (see
#' \code{\link{.discards_pretreatment_rows}}) comparable, in row count and CV
#' fold structure, to pipelines that do (\code{gene_level_fc}, which collapses
#' paired rows to one row per individual, or \code{rise_paired}/
#' \code{dearseq_mode = "paired"}, which model on post-treatment rows only).
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
