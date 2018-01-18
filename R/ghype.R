#' Title
#'
#' @param object
#' @param directed
#' @param selfloops
#' @param xi
#' @param omega
#' @param unbiased
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ghype <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, ...){
  UseMethod('ghype')
}

#' Title
#'
#' @param object
#' @param directed
#' @param selfloops
#' @param xi
#' @param omega
#' @param unbiased
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ghype.matrix <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, ...){
  if(is.null(xi)){
    xi=ComputeXi(object,directed,selfloops)
  }
  if(is.null(omega)){
    if(unbiased){
      omega <- matrix(1,nrow(object), ncol(object))
    } else{
      omega <- FitOmega(adj = object, xi = xi, directed = directed, selfloops = selfloops)
    }
  }

  model <- as.ghype(list('xi'= xi,
                         'omega' = omega,
                         'directed' = directed,
                         'selfloops' = selfloops))
  return(model)
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
as.ghype <- function(object, ...){
  UseMethod('as.ghype')
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
as.ghype.list <- function(object, ...){
  model <- list(
    'xi'= object$xi,
    'omega' = object$omega,
    'directed' = object$directed,
    'selfloops' = object$selfloops
  )
  class(model) <- 'ghype'
  return(model)
}
