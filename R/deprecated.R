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

#' @rdname link_significance
#' @export
#'
linkSignificance <- function(graph, model, under=FALSE, log.p=FALSE, binomial.approximation = FALSE, give_pvals = FALSE){
  warning('This function is deprecated and will be removed in future versions. Use `link_significance()`')
  link_significance(graph, model, under, log.p, binomial.approximation, give_pvals)
}

#' @rdname compute_xi
#' @export
#'
ComputeXi <- function(adj, directed, selfloops, regular = FALSE) {
  warning('This function is deprecated and will be removed in future versions. Use `compute_xi()`')
  compute_xi(adj, directed, selfloops, regular)
}