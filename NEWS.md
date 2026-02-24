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
