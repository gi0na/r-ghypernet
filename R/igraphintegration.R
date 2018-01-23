#' Title
#'
#' @param adjlist
#' @param directed
#' @param selfloops
#'
#' @return
#' @export
#'
#' @examples
CreateIgGraphs <- function(adjlist, directed, selfloops){
  if(directed)
    mode <- 'directed'
  if(!directed)
    mode <- 'undirected'

  lapply(X = adjlist, FUN = igraph::graph_from_adjacency_matrix, mode=mode, diag=selfloops)
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
ghype.igraph <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, ...){

  adj <- igraph::get.adjacency(graph = object, type = "both", sparse = FALSE)

  if(is.null(xi)){
    xi=ComputeXi(adj,directed,selfloops)
  }
  if(is.null(omega)){
    if(unbiased){
      omega <- matrix(1,nrow(adj), ncol(adj))
    } else{
      omega <- FitOmega(adj = adj, xi = xi, directed = directed, selfloops = selfloops)
    }
  }

  n <- nrow(adj)
  m <- sum(adj[mat2vec.ix(adj, directed, selfloops)])

  model <- as.ghype(list('adj' = adj,
                         'xi'= xi,
                         'omega' = omega,
                         'n' = n,
                         'm' = m,
                         'directed' = directed,
                         'selfloops' = selfloops))
  return(model)
}

#' Title
#'
#' @param graph
#' @param property
#' @param directed
#' @param selfloops
#' @param nsamples
#' @param xi
#' @param omega
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
BootstrapProperty <- function(graph, property, directed, selfloops, nsamples=1000, xi=NULL, omega=NULL, seed=NULL, ...){

  functionslist <- c(
    'page_rank',
    'page.rank',
    'page_rank_old',
    'page.rank.old'
  )

  if(directed){
    mode <- 'directed'
  } else{
    mode <- 'undirected'
  }

  if(!igraph::is.igraph(graph)){
    graph <- igraph::graph_from_adjacency_matrix(graph, mode=mode, diag = selfloops)
  }

  model <- ghype(object = graph, directed, selfloops, xi, omega)
  rsamples <- RandomGraph(nsamples, model, length(igraph::E(g)), seed=seed)
  gsamples <- CreateIgGraphs(adjlist = rsamples, directed = directed, selfloops = selfloops)
  if(as.character(substitute(property)) %in% functionslist){
    dproperty <- sapply(X = gsamples, FUN = function(graph, ...){match.fun(FUN = property)(graph, ...)$vector}, ...)
  } else{
    dproperty <- sapply(X = gsamples, FUN = property, ...)
  }
  return(dproperty)
}
