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
#' @param regular a boolean argument specifying whether to model the 'gnp' ensemble (no xi), defaults to FALSE.
#' @param ... additional arguments to be passed to the low level fitting functions.
#'
#' @return
#' ghype return an object of class "ghype".
#'
#' @export
#'
#'
ghype <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, regular=FALSE, ...){
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
#' @param regular a boolean argument specifying whether to model the 'gnp' ensemble (no xi), defaults to FALSE.
#' @param ... additional arguments to be passed to the low level fitting functions.
#'
#' @return
#' ghype return an object of class "ghype".
#'
#' @export
#'
#'
ghype.matrix <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, regular=FALSE, ...){

  df <- NULL

  if(is.null(xi)){
    xi=ComputeXi(object,directed,selfloops, regular = regular)
    df <- xiregular + (1-xiregular)*nrow(xi)*(1+directed)
  }

  if(is.null(omega)){
    if(unbiased){
      omega <- matrix(1,nrow(object), ncol(object))
    } else{
      omega <- FitOmega(adj = object, xi = xi, directed = directed, selfloops = selfloops)
      df <- df + sum(mat2vec.ix(omega,directed,selfloops)) - 1
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
#' @param regular a boolean argument specifying whether to model the 'gnp' ensemble (no xi), defaults to FALSE.
#' @param ... additional arguments to be passed to the low level fitting functions.
#'
#' @return
#' ghype return an object of class "ghype".
#'
#' @export
#'
#'
ghype.default <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, regular = FALSE, ...){

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

#' Fit the Soft-Configuration Model
#'
#' scm is wrapper for \link{ghype} that allows to specify a soft-configuration model.
#'
#' @param object either an adjacency matrix or an igraph graph
#' @param directed optional boolean, if not specified detected from object
#' @param selfloops optional boolean, if not specified detected from object
#' @param ... additional parameters
#'
#' @return ghype object
#' @export
#'
scm <- function(object, directed = NULL, selfloops = NULL, ...){

  specs <- check_specs(object)
  directed <- specs[1]
  selfloops <- specs[2]

  if(is.matrix(object)){
      if(!directed & !isSymmetric(object)){
        warning('Trying to compute undirected ensemble for asymmetric adjacency matrix.
              Adjacency matrix symmetrised as adj <- adj + t(adj)')
        object <- object + t(object)
      }
  }

  model <- ghype(object, directed=directed, selfloops=selfloops, unbiased = TRUE)
  model$df <- nrow(model$xi)*(1+directed)
  return(model)
}

#' Fit the gnm model
#'
#' regularm is wrapper for \link{ghype} that allows to specify a gnm regular model.
#' i.e. where all entries of the combinatorial matrix Xi are the same.
#'
#' @param object either an adjacency matrix or an igraph graph
#' @param directed optional boolean, if not specified detected from object
#' @param selfloops optional boolean, if not specified detected from object
#' @param ... additional parameters
#'
#' @return ghype object
#' @export
#'
regularm <- function(object, directed = NULL, selfloops = NULL, ...){

  specs <- check_specs(object)
  directed <- specs[1]
  selfloops <- specs[2]

  if(is.matrix(object)){
    if(!directed & !isSymmetric(object)){
      warning('Trying to compute undirected ensemble for asymmetric adjacency matrix.
              Adjacency matrix symmetrised as adj <- adj + t(adj)')
      object <- object + t(object)
    }
  }

  model <- ghype(object, directed=directed, selfloops=selfloops, unbiased = TRUE, regular = TRUE)
  model$df <- 1
  return(model)
}
