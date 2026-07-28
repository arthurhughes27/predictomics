# predictomics <img src="predictomics_hex.png" align="right" width="180" />

R package providing a standardised interface for prediction tasks in omics (e.g. transcriptomics) data. This package wraps the typical steps in predictive analyses with high-dimensional data (feature engineering, feature selection, model fitting, tuning, principled cross-validation) into a consistent, well-documented pipeline, so that different analysis choices can be swapped in and compared easily.

The core function is `predict_cv()`, a single cross-validated pipeline covering three stages, each swappable independently:

- **Feature engineering**: perform feature-level transformation (e.g. z-score, pre-post fold-change), geneset aggregation into pathway-level scores (mean, median, sum, first principal component, ssGSEA, GSVA).
- **Feature selection**: performs feature selection with different algorithms, including variance filtering, univariate correlation (Pearson/Spearman), cross-validated relative gain over a baseline, surrogate screening with the RISE algorithm, and `dearseq` differential expression filtering.
- **Model fitting**: Fits models according to a number of supported options `lm`, `glmnet`, `ridge`, `lasso`, `ranger` (random forest), and `svr`, each with inner-CV hyperparameter tuning via `caret`.

All three stages are re-fit within each outer CV fold by default to avoid data leakage. `predict_cv()` also supports treatment and covariate adjustment, k-fold or leave-one-out CV, and a covariates-only **baseline model** by simply passing `X = NULL`, giving a built-in null comparator for any analysis.

Beyond a single pipeline, `compare_pipelines()` benchmarks a set of alternative feature selection, engineering, or model choices against a baseline and a user-specified reference pipeline in one call, fitting everything sequentially and reporting RMSE/sRMSE/R²/Spearman correlation for each. 
