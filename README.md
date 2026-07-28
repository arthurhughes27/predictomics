# predictomics <img src="predictomics_hex.png" align="right" width="180" />

R package providing a standardised interface for prediction tasks in omics (e.g. transcriptomics) data. It wraps the usual mess of ad hoc scripts — filtering genes, aggregating them into pathway scores, fitting a model, tuning it, cross-validating it correctly — into one consistent, well-documented pipeline, so that different analysis choices can be swapped in and compared without rewriting boilerplate each time.

The core function is `predict_cv()`, a single cross-validated pipeline covering three stages, each swappable independently:

- **Feature engineering**: column transforms, geneset aggregation into pathway-level scores (mean, median, sum, first principal component, ssGSEA, GSVA), and paired pre/post-treatment fold-change collapsing (`gene_level_fc`).
- **Feature selection**: variance filtering, univariate correlation (Pearson/Spearman), cross-validated relative gain over a baseline, RISE surrogate screening, and `dearseq` differential expression filtering (gene- or geneset-level, classic or paired designs).
- **Model fitting**: `lm`, `glmnet`, `ridge`, `lasso`, `ranger` (random forest), and `svr`, each with inner-CV hyperparameter tuning via `caret`.

All three stages are re-fit within each outer CV fold by default to avoid data leakage (with an explicit, clearly-flagged `outside_cv` escape hatch for exploratory work). `predict_cv()` also supports treatment and covariate adjustment, k-fold or leave-one-out CV, and a covariates-only (or mean-only) **baseline model** by simply passing `X = NULL`, giving a built-in null comparator for any analysis.

Beyond a single pipeline, `compare_pipelines()` benchmarks a set of alternative filtration, engineering, or model choices against a baseline and a user-specified reference pipeline in one call, fitting everything sequentially and reporting RMSE/sRMSE/R²/Spearman correlation for each. Results come with `print()` and `plot()` methods styled consistently with the rest of the package, and `plot_selection_stability()` visualises how stable feature selection is across CV folds — useful for judging whether a filtering method is picking up signal or noise.
