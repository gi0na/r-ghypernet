#' Convert a list of adjacency matrices to a list of igraph graphs.
#'
#' @param adjlist a list of adjacency matrices
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#'
#' @return
#'
#' list of igraph graphs.
#'
#' @export
#'
CreateIgGraphs <- function(adjlist, directed, selfloops, weighted=FALSE){
  if(directed)
    mode <- 'directed'
  if(!directed)
    mode <- 'undirected'

  lapply(X = adjlist, FUN = igraph::graph_from_adjacency_matrix, mode=mode, diag=selfloops, weighted=weighted)
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
ghype.igraph <- function(object, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, ...){

  if(igraph::is_bipartite(object)){
    adj <- igraph::get.incidence(graph = object, sparse = FALSE)
  } else{
    adj <- igraph::get.adjacency(graph = object, type = "upper", sparse = FALSE)
    adj <- adj + t(adj)
  }

  if(is.null(xi)){
    xi=ComputeXi(adj,directed,selfloops)
  }

  df <- NULL

  if(is.null(omega)){
    df <- sum(mat2vec.ix(xi,directed,selfloops))
    if(unbiased){
      omega <- matrix(1,nrow(adj), ncol(adj))
    } else{
      omega <- FitOmega(adj = adj, xi = xi, directed = directed, selfloops = selfloops)
      df <- df + sum(mat2vec.ix(omega,directed,selfloops))
    }
  }

  if(nrow(adj)==ncol(adj)){
    n <- nrow(adj)
  } else{
    n <- nrow(adj)+ncol(adj)
  }

  m <- sum(adj[mat2vec.ix(adj, directed, selfloops)])

  model <- as.ghype(list(call = match.call(),
                         'adj' = adj,
                         'xi'= xi,
                         'omega' = omega,
                         'n' = n,
                         'm' = m,
                         'directed' = directed,
                         'selfloops' = selfloops,
                         'df' = df ))
  return(model)
}

#' TODO
#'
#' .....
#'
#' @param graph adj
#' @param property function
#' @param directed ..
#' @param selfloops ..
#' @param nsamples ..
#' @param xi ..
#' @param omega ..
#' @param ... ...
#' @param m ...
#' @param seed ..
#'
#' @return vector
#' @export
#'
BootstrapProperty <- function(graph, property, directed, selfloops, nsamples=1000, xi=NULL, omega=NULL, model=NULL, m=NULL, seed=NULL, ...){

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
    if(nrow(graph)!=ncol(graph)){
      # default to out direction for bipartite graphs
      graph <- igraph::graph_from_incidence_matrix(graph, directed = directed, mode='out')
    } else{
      graph <- igraph::graph_from_adjacency_matrix(graph, mode=mode, diag = selfloops)
    }
  }
  if(is.null(m))
    m <- length(igraph::E(graph))

  if(is.null(model))
    model <- ghype(object = graph, directed, selfloops, xi, omega)

  rsamples <- RandomGraph(nsamples, model, m, seed=seed)
  gsamples <- CreateIgGraphs(adjlist = rsamples, directed = directed, selfloops = selfloops)
  if(as.character(substitute(property)) %in% functionslist){
    dproperty <- sapply(X = gsamples, FUN = function(graph, directed, ...){match.fun(FUN = property)(graph, directed=directed, ...)$vector}, directed=directed, ...)
  } else{
    dproperty <- sapply(X = gsamples, FUN = property, ...)
  }
  return(dproperty)
}
