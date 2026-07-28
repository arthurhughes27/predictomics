# predictomics <img src="predictomics_hex.png" align="right" width="180" />

R package providing a standardised interface for prediction tasks in omics data.

Built around `predict_cv()`, a single cross-validated pipeline function covering feature engineering (fold-change transforms, geneset aggregation via mean/median/sum/pc1/ssGSEA/GSVA), feature selection (variance, correlation, relative gain, RISE, dearseq differential expression), and model fitting (`lm`, `glmnet`, `ridge`, `lasso`, `ranger`, `svr`) with inner-CV hyperparameter tuning. Supports treatment/covariate adjustment, paired pre/post-treatment designs, and a covariates-only baseline model.

`compare_pipelines()` benchmarks alternative filtration, engineering, or model choices against a baseline and reference pipeline in one call. Results can be inspected via `print()`/`plot()` methods, and `plot_selection_stability()` visualises feature selection stability across CV folds.
