# SAGA2 2.0.1

## Bug fixes

* `EvaluateModel()` now computes the Bonferroni-corrected alpha using the
  total number of models evaluated, as documented. It previously divided
  0.05 by the number of elements in the `genarch` object (always 7), so
  reported alpha levels were far too permissive.
* `LCA()`: the `best.eqns.w` element is now correctly aligned with
  `best.models` (both in decreasing weight order). Previously the weights
  were returned in equation-index order, so pairing weights with models
  gave wrong values.
* Data validation now catches any invalid sire/dam entry (out of range,
  non-integer, or NA) and names the offending rows. Previously a single
  bad value passed validation and the affected cohorts were silently
  dropped from all model fits.
* Fixes applied during validation (the `X0`/`Z0` spelling of the sex
  chromosome system, all-`NA` environmental columns, and lowercase `se`
  column names) now propagate to the analysis instead of being lost;
  `LCA(data, SCS = "X0")` and all-NA `enviro` columns no longer error.
* Internal duplicate-column detection in the C-matrix no longer fails on
  minimal (two CGE column) C-matrices.
* `plotObserved()` now fits the additive expectation line with weights
  `SE^-2` (previously raw `SE`, which gave the most weight to the least
  precise cohorts) and fits on the unjittered predictor, so the line no
  longer moves between calls.
* `EvaluateModel()` restores graphical parameters on exit and errors
  informatively when the `genarch` object was created with
  `ret.all = FALSE`.
* The message reporting how many models were removed due to high
  covariance no longer over-counts by one.
* Aliased coefficients are now dropped by their `NA` status rather than
  by position, so user-supplied C-matrices are handled robustly; a user
  C-matrix must (and is now checked to) contain an all-1s mean column.

## Behavior changes

* **AICc is now computed with the standard small-sample formula**
  `AICc = AIC + 2K(K+1)/(n - K - 1)` with `K` equal to the number of
  estimated parameters including the residual variance. The previous
  formula used `2k(k+1)/(n - k)` with `k` one less. Numeric results
  (delta AICc, Akaike weights, model-averaged estimates, unconditional
  SEs, and variable importance) therefore differ from version 2.0.0 and
  model rankings can change for small numbers of cohorts.
* `$estimates` and `$varimp` are now numeric matrices (previously
  character). `$varimp` is a one-column matrix with CGE names as row
  names.
* The `Estiamtes` element of `EvaluateModel()` output has been renamed
  `Estimates`; `Estiamtes` is retained as a deprecated duplicate alias
  and will be removed in a future release.

## Performance

* The exhaustive model-scoring loop now uses `lm.wfit()` instead of
  `glm()` (identical weighted least squares results, verified to
  reproduce the glm AICc values exactly), computes weights once instead
  of once per model, and vectorizes variable importance with `rowsum()`.
  On the shipped BanInf data `LCA()` is roughly 3x faster; the advantage
  grows with model space size. Fitted `glm` objects are still returned
  for the confidence set and, when `ret.all = TRUE`, for all models.

# SAGA2 2.0.0

* Initial release of SAGA2 with automatic C-matrix construction from user-supplied breeding designs.
* Main analysis function `LCA()` replaces the previous `AnalyzeCrossesMM()` interface.
* Added support for environmental variation and gene-by-environment interactions via the `env` argument.
* Added observed phenotypic parental effects as an alternative to calculated maternal effects via `parental = "obs"`.
* Added support for XY, XO, ZW, ZO, and NSC sex chromosome systems.
* Added `plotObserved()` for traditional line cross analysis plots.
* Added `VisModelSpace()` for visualizing Akaike weight distributions across model space.
* Added `EvaluateModel()` for examining individual model results.
* S3 `plot()` method for `genarch` objects provides publication-quality figures with variable importance color coding.
