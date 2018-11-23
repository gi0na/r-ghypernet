updateModel <- function(model, adj){

  callname <- as.character(model$call[1])

  xi <- NULL

  newcall <- NULL
  if(length(grep('function', callname))>0){
    if(length(grep('block', callname))>0 | length(grep('labels', callname))>0){
      callname <- 'bccm'
    } else{
      callname <- 'ghype'
    }
  } else{
    fixXi <- length(grep('xi', deparse(model$call, width.cutoff = 500)))>0
    if(fixXi)
      xi <- model$xi
  }
  if(length(grep('ghype', callname))>0){
    callname <- 'ghype'
    newcall <- call(name = callname, object=adj, directed=model$directed, selfloops=model$selfloops, xi=xi, unbiased=all(model$omega==1))
  } else{
    if(length(grep('bccm', callname))>0){
      callname <- 'bccm'
      newcall <- call(name = callname, adj=adj, labels=model$labels, directed=model$directed, selfloops=model$selfloops, xi=xi)
    }
    if(length(grep('nrm', callname))>0){
      callname <- 'ghype'
      newcall <- call(name = callname, object=adj, directed=model$directed, selfloops=model$selfloops, xi=xi, omega=model$omega)
    }
  }

  return(newcall)
}

#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
AIC.ghype <- function(object, ...){
  2*object$df - 2*object$loglikelihood
}

#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
BIC.ghype <- function(object, ...){
  log(object$m)*object$df - 2*object$loglikelihood
}
