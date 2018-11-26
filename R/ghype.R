#' Fitting gHypEG models
#'
#' ghype is used to fit gHypEG models when the propensity matrix is known.
#' It can be used to estimate a null model (soft configuration model), or
#' the benchmark 'full-model', where the propensity matrix is fitted such
#' that the expected graph from the fitted model is the one passed to the
#' function.
#'
#' @param object either an adjacency matrix or an igraph graph.
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param xi an optional matrix defining the combinatorial matrix of the model.
#' @param omega an optional matrix defining the propensity matrix of the model.
#' @param unbiased a boolean argument specifying whether to model the hypergeometric ensemble (no propensity), defaults to FALSE.
#' @param ... additional arguments to be passed to the low level fitting functions.
#'
#' @return
#' ghype return an object of class "ghype".
#'
#' @export
#'
#'
ghype <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, ...){
  UseMethod('ghype')
}


#' Fitting gHypEG models
#'
#' ghype is used to fit gHypEG models when the propensity matrix is known.
#' It can be used to estimate a null model (soft configuration model), or
#' the benchmark 'full-model', where the propensity matrix is fitted such
#' that the expected graph from the fitted model is the one passed to the
#' function.
#'
#' @param object either an adjacency matrix or an igraph graph.
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param xi an optional matrix defining the combinatorial matrix of the model.
#' @param omega an optional matrix defining the propensity matrix of the model.
#' @param unbiased a boolean argument specifying whether to model the hypergeometric ensemble (no propensity), defaults to FALSE.
#' @param ... additional arguments to be passed to the low level fitting functions.
#'
#' @return
#' ghype return an object of class "ghype".
#'
#' @export
#'
#'
ghype.matrix <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, ...){
  if(is.null(xi)){
    xi=ComputeXi(object,directed,selfloops)
  }

  df <- NULL

  if(is.null(omega)){
    df <- sum(mat2vec.ix(xi,directed,selfloops))
    if(unbiased){
      omega <- matrix(1,nrow(object), ncol(object))
    } else{
      omega <- FitOmega(adj = object, xi = xi, directed = directed, selfloops = selfloops)
      df <- df + sum(mat2vec.ix(omega,directed,selfloops))
    }
  }

  if(nrow(object)==ncol(object)){
    n <- nrow(object)
  } else{
    n <- c(nrow(object)+ncol(object),nrow(object),ncol(object))
  }

  m <- sum(object[mat2vec.ix(object, directed, selfloops)])

  model <- as.ghype(list(call = match.call(),
                         'adj' = object,
                         'xi'= xi,
                         'omega' = omega,
                         'n' = n,
                         'm' = m,
                         'directed' = directed,
                         'selfloops' = selfloops,
                         'df' = df))
  return(model)
}


#' Fitting gHypEG models
#'
#' ghype is used to fit gHypEG models when the propensity matrix is known.
#' It can be used to estimate a null model (soft configuration model), or
#' the benchmark 'full-model', where the propensity matrix is fitted such
#' that the expected graph from the fitted model is the one passed to the
#' function.
#'
#' @param object either an adjacency matrix or an igraph graph.
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param xi an optional matrix defining the combinatorial matrix of the model.
#' @param omega an optional matrix defining the propensity matrix of the model.
#' @param unbiased a boolean argument specifying whether to model the hypergeometric ensemble (no propensity), defaults to FALSE.
#' @param ... additional arguments to be passed to the low level fitting functions.
#'
#' @return
#' ghype return an object of class "ghype".
#'
#' @export
#'
#'
ghype.default <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, ...){

  if(is.null(omega)){
    if(unbiased){
      omega <- matrix(1,nrow(object), ncol(object))
    }
  }

  n <- nrow(xi)
  m <- sqrt(sum(xi))

  model <- as.ghype(list(call = match.call(),
                         'adj' = object,
                         'xi'= xi,
                         'omega' = omega,
                         'n' = n,
                         'm' = m,
                         'directed' = directed,
                         'selfloops' = selfloops,
                         'df' = NULL))
  return(model)
}


#' Map list to ghype object
#'
#' Manually map a list to a ghype object
#'
#' @param object list object to map to ghype.
#' @param ... additional arguments to be passed to the low level functions.
#'
#' @return
#' an object of class "ghype"
#'
#' @export
#'
as.ghype <- function(object, ...){
  UseMethod('as.ghype')
}


#' Map list to ghype object
#'
#' Manually map a list to a ghype object
#'
#' @param object list object to map to ghype.
#' @param ... additional arguments to be passed to the low level functions.
#'
#' @return
#' an object of class "ghype"
#'
#' @export
#'
as.ghype.list <- function(object, ...){
  model <- list(
    call = object$call,
    'adj' = object$adj,
    'n' = object$n,
    'm' = object$m,
    'xi'= object$xi,
    'omega' = object$omega,
    'directed' = object$directed,
    'selfloops' = object$selfloops,
    'loglikelihood' = object$loglikelihood,
    'df' = object$df
  )
  if(is.null(model$loglikelihood) & !is.null(model$adj)){
    model$loglikelihood <- logl(adj=model$adj, xi=model$xi,
                                omega=model$omega, directed=model$directed,
                                selfloops=model$selfloops)
  }
  class(model) <- 'ghype'
  return(model)
}
