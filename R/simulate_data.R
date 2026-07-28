# =============================================================================
# simulate_data.R
# Simulated omics data generator for examples, tests, and the package
# tutorial vignette.
# =============================================================================


# -----------------------------------------------------------------------------
#' Simulate gene expression data for predictomics examples
#'
#' @description
#' Generates a simulated gene-expression-like dataset with a known signal,
#' for use in examples, tests, and the package tutorial vignette. A handful
#' of genesets are "signal" genesets whose mean expression drives the
#' response \code{Y}; the rest are pure noise. Also simulates a binary
#' \code{treatment}, covariates (\code{age}, \code{sex}), and, optionally, a
#' paired pre/post-treatment design suitable for
#' \code{engineering_params$gene_level_fc}, \code{selection_params$dearseq_mode
#' = "paired"}, or \code{selection_params$rise_paired}.
#'
#' @details
#' Genes are grouped into \code{n_genesets} non-overlapping genesets of
#' \code{geneset_size} genes each, named \code{"geneset1"}, \code{"geneset2"},
#' etc. The first \code{n_signal_genesets} are signal genesets: their member
#' genes are simulated with a shared per-sample latent factor added, so that
#' the geneset mean is informative for \code{Y}. Remaining genes are pure
#' Gaussian noise, uncorrelated with \code{Y}.
#'
#' For \code{design = "cross_sectional"} (default), \code{n} independent
#' samples are simulated, \code{Y} is a linear combination of the signal
#' geneset latent factors (plus a \code{treatment} effect and noise), and a
#' binary \code{treatment} and covariates (\code{age}, \code{sex}) are
#' generated independently of the signal.
#'
#' For \code{design = "paired"}, \code{n} individuals are each simulated at
#' \code{timepoint = 0} (pre) and \code{timepoint = 1} (post). The
#' post-treatment change in the signal genesets' latent factors drives
#' \code{Y} (measured post-treatment only), mimicking a pre/post intervention
#' study; \code{treatment} encodes two study arms, independent of the
#' pre/post structure.
#'
#' @param n Integer. Number of samples (\code{design = "cross_sectional"}) or
#'   individuals (\code{design = "paired"}, giving \code{2 * n} rows).
#'   Defaults to \code{60}.
#' @param p Integer. Total number of genes. Must be a multiple of
#'   \code{geneset_size}. Defaults to \code{200}.
#' @param n_genesets Integer. Number of genesets to partition the \code{p}
#'   genes into. Defaults to \code{10}.
#' @param geneset_size Integer. Number of genes per geneset. Defaults to
#'   \code{20}.
#' @param n_signal_genesets Integer. Number of genesets (from the first
#'   \code{n_genesets}) whose latent factor drives \code{Y}. Defaults to
#'   \code{2}.
#' @param design Character string. One of \code{"cross_sectional"} (default)
#'   or \code{"paired"}.
#' @param treatment_effect Numeric. Size of the \code{treatment} effect added
#'   to \code{Y} (cross-sectional design only). Defaults to \code{1}.
#' @param noise_sd Numeric. Standard deviation of the noise added to \code{Y}.
#'   Defaults to \code{1}.
#' @param seed Integer. Random seed for reproducibility. Defaults to \code{1}.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{X}}{Numeric matrix of gene expression (samples x genes).}
#'     \item{\code{Y}}{Numeric response vector.}
#'     \item{\code{treatment}}{Binary numeric vector (0/1).}
#'     \item{\code{covariates}}{Data frame with columns \code{age} (numeric)
#'       and \code{sex} (factor).}
#'     \item{\code{genesets}}{Named list of character vectors of gene names.}
#'     \item{\code{signal_genesets}}{Character vector of the names of the
#'       genesets that actually drive \code{Y} (useful for checking that
#'       feature selection recovers the right signal).}
#'     \item{\code{individual_id}, \code{timepoint}}{Only present when
#'       \code{design = "paired"}.}
#'   }
#'
#' @seealso \code{\link{predict_cv}}, \code{\link{compare_pipelines}}
#'
#' @examples
#' d <- simulate_predictomics_data(n = 40, p = 100)
#' dim(d$X)
#' d$signal_genesets
#'
#' @export
# -----------------------------------------------------------------------------
simulate_predictomics_data <- function(n                 = 60,
                                       p                 = 200,
                                       n_genesets        = 10,
                                       geneset_size      = 20,
                                       n_signal_genesets = 2,
                                       design            = "cross_sectional",
                                       treatment_effect  = 1,
                                       noise_sd          = 1,
                                       seed              = 1) {

  if (!design %in% c("cross_sectional", "paired"))
    stop("[predictomics] design must be 'cross_sectional' or 'paired'.",
         call. = FALSE)

  if (p != n_genesets * geneset_size)
    stop("[predictomics] p must equal n_genesets * geneset_size.",
         call. = FALSE)

  set.seed(seed)

  gene_names <- paste0("gene", seq_len(p))
  genesets   <- split(gene_names, rep(paste0("geneset", seq_len(n_genesets)),
                                      each = geneset_size))
  signal_genesets <- names(genesets)[seq_len(n_signal_genesets)]

  if (design == "cross_sectional") {

    # One latent factor per signal geneset, shared by its member genes
    latent <- matrix(rnorm(n * n_signal_genesets), nrow = n, ncol = n_signal_genesets)
    X <- matrix(rnorm(n * p), nrow = n, ncol = p, dimnames = list(NULL, gene_names))
    for (i in seq_len(n_signal_genesets)) {
      genes <- genesets[[signal_genesets[i]]]
      X[, genes] <- X[, genes] + latent[, i]
    }

    treatment  <- rbinom(n, 1, 0.5)
    covariates <- data.frame(
      age = round(rnorm(n, mean = 45, sd = 12)),
      sex = factor(sample(c("F", "M"), n, replace = TRUE))
    )

    Y <- rowSums(latent) + treatment_effect * treatment + rnorm(n, sd = noise_sd)

    rownames(X) <- paste0("sample", seq_len(n))

    list(
      X               = X,
      Y               = Y,
      treatment       = treatment,
      covariates      = covariates,
      genesets        = genesets,
      signal_genesets = signal_genesets
    )

  } else {

    n_total       <- n * 2
    individual_id <- rep(seq_len(n), each = 2)
    timepoint     <- rep(c(0, 1), times = n)

    # Individual-level baseline expression, plus a post-treatment shift in
    # the signal genesets that varies across individuals (the "fold change")
    baseline_shift <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
    fc             <- matrix(rnorm(n * n_signal_genesets), nrow = n, ncol = n_signal_genesets)

    X <- matrix(rnorm(n_total * p), nrow = n_total, ncol = p,
               dimnames = list(NULL, gene_names))
    X <- X + baseline_shift[individual_id, ]

    for (i in seq_len(n_signal_genesets)) {
      genes    <- genesets[[signal_genesets[i]]]
      post_idx <- which(timepoint == 1)
      X[post_idx, genes] <- X[post_idx, genes] + fc[, i]
    }

    treatment  <- rep(rbinom(n, 1, 0.5), each = 2)
    covariates <- data.frame(
      age = rep(round(rnorm(n, mean = 45, sd = 12)), each = 2),
      sex = factor(rep(sample(c("F", "M"), n, replace = TRUE), each = 2))
    )

    # Y (measured post-treatment) driven by the size of the post-pre shift.
    # predict_cv()'s gene_level_fc mode only ever uses the post-treatment Y
    # value, but Y must still be fully non-NA at every row, so the
    # (unused) pre-treatment entry is just set equal to the post value.
    Y_post <- rowSums(fc) + rnorm(n, sd = noise_sd)
    Y                 <- numeric(n_total)
    Y[timepoint == 1] <- Y_post
    Y[timepoint == 0] <- Y_post

    rownames(X) <- paste0("sample", seq_len(n_total))

    list(
      X               = X,
      Y               = Y,
      treatment       = treatment,
      covariates      = covariates,
      individual_id   = individual_id,
      timepoint       = timepoint,
      genesets        = genesets,
      signal_genesets = signal_genesets
    )
  }
}
