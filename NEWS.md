# ghypernet 1.1.0

- function renaming:
* `linkSignificance()` --> `link_significance()`
* `ComputeXi()` --> `compute_xi()`

- new features
* updated documentation for endogenous statistics functions
* introduced a new function `get_zero_dummy()` that generates the dummy variable to be used together with endogenous statistics functions.
* added testing with test_that
* `check_specs()`: recognize bipartite graphs
* addressed issue #15: bccm now names coefficients according to label names. E.g.: coefficient for the propensities between block 'red' and 'blue' are now named 'red<->blue'.

- bug fixes:
* `scm()`: fixed df for bipartite graphs
* `linkSignificance()`: fixed issue for bipartite graphs
* `bccm()`: fixed issue when fitting to bipartite graphs
* `bccm()`: fixed issue arising when passing a list of identical labels (all vertices in one block)

# ghypernet 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Fixed a bug in `check_specs()` that threw an error when the function was called without additional parameters.
* Modified `bccm()` to improve its speed when passing a large number of labels.
* Renamed `nrmSelection()`, `createPredictors`, and `nrmChoose()` to `nrm_selection()`, `create_predictors`, and `nrm_choose()` respectively, to standardise functions' names.
* Updated `nrm_selection()` and `nrm_choose()` to use 'pbmclapply' and show a progress bar.
* Added compatibility to 'texreg' to visualise the output of nrm as a table of parameters.

