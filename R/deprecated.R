#' @rdname nrm_selection
#' @export
nrmSelection <- function(adj, predictors, 
                          directed, selfloops, pval = 0.05, 
                          xi = NULL, init = NULL, ncores = NULL, 
                          ...){
  warning('This function is deprecated and will be removed in future versions. Use `nrm_selection()`')
  nrm_selection(adj, predictors, 
                directed, selfloops, pval, 
                xi, init, ncores, 
                ...)
}

#' @rdname nrm_choose
#' @export
nrmChoose <- function(adj, w.list, 
                       xi = NULL, directed, selfloops, 
                       pval = 0.05, init = NULL, ncores = NULL) {
  warning('This function is deprecated and will be removed in future versions. Use `nrm_choose()`')
  nrm_choose(adj, w.list, 
                xi, directed, selfloops, 
                pval, init, ncores)
}

#' @rdname create_predictors
#' @export
createPredictors <- function(predictors, 
                              ...) {
  warning('This function is deprecated and will be removed in future versions. Use `create_predictors()`')
  create_predictors(predictors,...)
}