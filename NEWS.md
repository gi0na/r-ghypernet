# ghypernet 1.0.1.0921

- function renaming:
* `linkSignificance()` --> `link_significance()`

- new features
* `check_specs()`: recognize bipartite graphs

- bug fixes:
* `scm()`: fixed df for bipartite graphs
* `linkSignificance()`: fixed issue for bipartite graphs
* `bccm()`: fixed issue for bipartite graphs

# ghypernet 1.0.1.0910

* Fixed bug in `phi_element_wallenius()`: BiasedUrn throws an error when there is a highly improbable observation. Substitutes CDF values (error) with 1 (good approx).

# ghypernet 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Fixed a bug in `check_specs()` that threw an error when the function was called without additional parameters.
* Modified `bccm()` to improve its speed when passing a large number of labels.
* Renamed `nrmSelection()`, `createPredictors`, and `nrmChoose()` to `nrm_selection()`, `create_predictors`, and `nrm_choose()` respectively, to standardise functions' names.
* Updated `nrm_selection()` and `nrm_choose()` to use 'pbmclapply' and show a progress bar.
* Added compatibility to 'texreg' to visualise the output of nrm as a table of parameters.

