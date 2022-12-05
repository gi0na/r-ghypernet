# ghypernet 1.1.2

- bug fixes:
* `bccm()`: fixed bug that returned wrongly filled omegaBlock matrix. The bug should not have affected the model fit itself.

# ghypernet 1.1.1

- new features:
* extended testing scripts of package

- bug fixes:
* `logl()` and `rghype()`: fixed bug that treated as a hypergeometric model a ghype where all odds are equal to 0. Throws error instead.
* `compute_xi()`: Fixed a bug that introduced spurious stubs for nodes with zero degrees when creating a Xi matrix without selfloops
* `link_significance()`: fixed bug when under=TRUE, updated implementation to improve performance
* `ghype()`: fixed computation of degrees of freedom for full model where xi has zero values

# ghypernet 1.1.0.1

- bug fixes:
* `link_significance()`: When calling with under = FALSE and pval = FALSE was returning phi(0) = 1 instead of Pr(x>0)

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

