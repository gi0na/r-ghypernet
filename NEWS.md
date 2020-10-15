# ghypernet 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Fixed a bug in `check_specs()` that threw an error when the function was called without additional parameters.
* Modified `bccm()` to improve its speed when passing a large number of labels.
* Renamed `nrmSelection()`, `createPredictors`, and `nrmChoose()` to `nrm_selection()`, `create_predictors`, and `nrm_choose()` respectively, to standardise functions' names.
* Updated `nrm_selection()` and `nrm_choose()` to use 'pbmclapply' and show a progress bar.
* Added compatibility to 'texreg' to visualise the output of nrm as a table of parameters.

