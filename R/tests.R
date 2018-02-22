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
llratiotest <- function(nullmodel, altmodel, df){
  llratio <- loglratio(nullmodel,altmodel)
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
  n <- full$n
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
