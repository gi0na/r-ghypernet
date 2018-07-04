#' Title
#'
#' @param nullmodel
#' @param altmodel
#' @param df
#'
#' @return
#' @export
#'
#' @examples
llratiotest <- function(nullmodel, altmodel, df=NULL){
  llratio <- loglratio(nullmodel,altmodel)
  if(is.null(df)){
    df <- altmodel$df-nullmodel$df
  }
  return(pchisq(q = -2*llratio, df = df, lower.tail = FALSE))
}

#' Title
#'
#' @param graph
#' @param directed
#' @param selfloops
#'
#' @return
#' @export
#'
#' @examples
isNetwork <- function(graph, directed, selfloops){
  full <- ghype(graph, directed, selfloops)
  null <- ghype(graph, directed, selfloops, unbiased = TRUE)
  n <- full$n[1]
  df <- n*(n-!selfloops)/(1+!directed)
  if(igraph::is.igraph(graph)){
    if(igraph::is.bipartite(graph))
      df <- sum(V(graph)$type)*sum(!V(graph)$type)
  }
  if(is.matrix(graph)){
    if(nrow(graph)!=ncol(graph)){
      df <- nrow(graph)*ncol(graph)
    }
  }
  return(llratiotest(null,full,df))
}

#' Title
#'
#' @param graph
#' @param directed
#' @param selfloops
#'
#' @return
#' @export
#'
#' @examples
linkSignificance <- function(graph, model){
  ## TODO: assume graph is adjacency for now. Extend with method for graph and for matrices

  directed <- model$directed
  selfloops <- model$selfloops

  # get relevant indices
  idx <- mat2vec.ix(graph, directed, selfloops)

  # compute parameters for marginal distributions
  xibar <- sum(model$xi[idx])-model$xi[idx]
  omegabar <- (sum(model$xi[idx]*model$omega[idx])-model$xi[idx]*model$omega[idx])/xibar

  # compute vector of probabilities using Wallenius univariate distribution
  probvec <- Vectorize(FUN = function(id){
    BiasedUrn::pWNCHypergeo(x = graph[idx][id],m1 = model$xi[idx][id],m2 = xibar[id], n = sum(graph[idx]), odds = model$omega[idx][id]/omegabar[id])
    }, vectorize.args = 'id')(1:sum(idx))

  # return matrix of significance for each entry of original adjacency
  return(vec2mat(probvec,directed,selfloops,nrow(graph)))
}
