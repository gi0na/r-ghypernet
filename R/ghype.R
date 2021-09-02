#' Fitting gHypEG models
#'
#' ghype is used to fit gHypEG models when the propensity matrix is known.
#' It can be used to estimate a null model (soft configuration model), or
#' the benchmark 'full-model', where the propensity matrix is fitted such
#' that the expected graph from the fitted model is the one passed to the
#' function.
#'
#' @param graph either an adjacency matrix or an igraph graph.
#' @param directed a boolean argument specifying whether graph is directed or not.
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
#' @examples
#' data("adj_karate")
#' fullmodel <- ghype(graph = adj_karate, directed = FALSE, selfloops = FALSE, unbiased = FALSE)
#'
ghype <- function(graph, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, regular=FALSE, ...){
  UseMethod('ghype')
}


#' @describeIn ghype Fitting ghype models from an adjacency matrix
#'
#' @export
#'
#'
ghype.matrix <- function(graph, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, regular=FALSE, ...){

  df <- NULL

  if(is.null(xi)){
    xi=ComputeXi(graph,directed,selfloops, regular = regular)
    df <- regular + (1-regular)*nrow(xi)*(1+directed)
  }

  if(is.null(omega)){
    if(unbiased){
      omega <- matrix(1,nrow(graph), ncol(graph))
    } else{
      omega <- FitOmega(adj = graph, xi = xi, directed = directed, selfloops = selfloops)
      df <- df + sum(mat2vec.ix(omega,directed,selfloops)) - 1
    }
  }

  if(nrow(graph)==ncol(graph)){
    n <- nrow(graph)
  } else{
    n <- c(nrow(graph)+ncol(graph),nrow(graph),ncol(graph))
  }

  m <- sum(graph[mat2vec.ix(graph, directed, selfloops)])

  model <- as.ghype(list(call = match.call(),
                         'adj' = graph,
                         'xi'= xi,
                         'omega' = omega,
                         'n' = n,
                         'm' = m,
                         'directed' = directed,
                         'selfloops' = selfloops,
                         'regular' = regular,
                         'unbiased' = unbiased,
                         'df' = df), ...)
  return(model)
}


#' @describeIn ghype Generating a ghype model from given xi and omega
#'
#' @export
#'
#'
ghype.default <- function(graph, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, regular = FALSE, ...){

  if(is.null(omega) & is.matrix(graph) & unbiased){
    omega <- matrix(1,nrow(graph), ncol(graph))
  }

  n <- nrow(xi)
  m <- sqrt(sum(xi))

  model <- as.ghype(list(call = match.call(),
                         'adj' = graph,
                         'xi'= xi,
                         'omega' = omega,
                         'n' = n,
                         'm' = m,
                         'directed' = directed,
                         'selfloops' = selfloops,
                         'regular' = regular,
                         'unbiased' = unbiased,
                         'df' = NULL))
  return(model)
}


#' Map list to ghype object
#'
#' Manually map a list to a ghype object
#'
#' @param object list object to map to ghype.
#' @param ... additional arguments to be passed to logl function.
#'
#' @return
#' an object of class "ghype"
#'
#' @export
#' 
#' @examples
#' ll <- list(call = NULL, 'adj' = NULL, 'xi'= matrix(36,4,4), 'omega' = matrix(1,4,4), 
#'      'n' = 4, 'm' = 12, 'directed' = TRUE, 'selfloops' = TRUE,
#'      'regular' = TRUE, 'unbiased' = TRUE, 'df' = 1)
#' model <- as.ghype(ll)
#'
as.ghype <- function(object, ...){
  UseMethod('as.ghype')
}


#' @describeIn as.ghype Map list to ghype
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
    'regular' = object$regular,
    'unbiased' = object$unbiased,
    'df' = object$df
  )
  if(is.null(model$loglikelihood) & !is.null(model$adj) & !is.null(model$xi) & !is.null(model$omega)){
    model$loglikelihood <- logl(object=model$adj, xi=model$xi,
                                omega=model$omega, directed=model$directed,
                                selfloops=model$selfloops, ...)
  }
  class(model) <- 'ghype'
  return(model)
}

#' Fit the Soft-Configuration Model
#'
#' scm is wrapper for \link{ghype} that allows to specify a soft-configuration model.
#'
#' @param graph either an adjacency matrix or an igraph graph
#' @param directed optional boolean, if not specified detected from graph
#' @param selfloops optional boolean, if not specified detected from graph
#' @param ... additional parameters passed to the ghype function
#'
#' @return ghype object
#' @export
#' 
#' @examples
#' data("adj_karate")
#' confmodel <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#'
scm <- function(graph, directed = NULL, selfloops = NULL, ...){

  if(is.null(directed) | is.null(selfloops)){
    specs <- check_specs(graph)
    if(is.null(directed)) directed <- specs[1]
    if(is.null(selfloops)) selfloops <- specs[2]
  }

  if(is.matrix(graph)){
      if(!directed & !isSymmetric(graph)){
        warning('Trying to compute undirected ensemble for asymmetric adjacency matrix.
              Adjacency matrix symmetrised as adj <- adj + t(adj)')
        graph <- graph + t(graph)
      }
  }

  model <- ghype(graph, directed=directed, selfloops=selfloops, unbiased = TRUE, regular = FALSE, ...)
  model$df <- nrow(model$xi) * (directed) + ncol(model$xi)
  return(model)
}

#' Fit the gnm model
#'
#' regularm is wrapper for \link{ghype} that allows to specify a gnm regular model.
#' i.e. where all entries of the combinatorial matrix Xi are the same.
#'
#' @param graph either an adjacency matrix or an igraph graph
#' @param directed optional boolean, if not specified detected from graph
#' @param selfloops optional boolean, if not specified detected from graph
#' @param ... additional parameters passed to the ghype function
#'
#' @return ghype object
#' @export
#' 
#' @examples
#' data("adj_karate")
#' regularmodel <- regularm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#'
regularm <- function(graph, directed = NULL, selfloops = NULL, ...){

  if(is.null(directed) | is.null(selfloops)){
    specs <- check_specs(graph)
    if(is.null(directed)) directed <- specs[1]
    if(is.null(selfloops)) selfloops <- specs[2]
  }

  if(is.matrix(graph)){
    if(!directed & !isSymmetric(graph)){
      warning('Trying to compute undirected ensemble for asymmetric adjacency matrix.
              Adjacency matrix symmetrised as adj <- adj + t(adj)')
      graph <- graph + t(graph)
    }
  }

  model <- ghype(graph, directed=directed, selfloops=selfloops, unbiased = TRUE, regular = TRUE, ...)
  model$df <- 1
  return(model)
}
