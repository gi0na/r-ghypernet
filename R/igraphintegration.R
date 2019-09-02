# auxiliary function for to extract igraph properties
check_specs.igraph <- function(object, ...){
  if(requireNamespace("igraph", quietly = TRUE) && igraph::is.igraph(object)){
    if(is.null(directed)){
      if(igraph::is.directed(object)){
        directed <- FALSE
      } else{
        directed <- TRUE
      }
    }

    if(is.null(selfloops)){
      if(igraph::is.simple(igraph::simplify(object, remove.multiple = TRUE, remove.loops = FALSE))){
        selfloops <- FALSE
      } else{
        selfloops <- TRUE
      }
    }
  }
  return(c('directed'=directed, 'selfloops'=selfloops))
}

#' Convert a list of adjacency matrices to a list of igraph graphs.
#'
#' @param adjlist a list of adjacency matrices
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param weighted boolean, generate weighted graphs?
#'
#' @return
#'
#' list of igraph graphs.
#'
#' @export
#' 
#' @examples
#' data('adj_karate')
#' adj_list <- list(adj_karate)
#' glist <- CreateIgGraphs(adj_list, FALSE, FALSE)
#'
CreateIgGraphs <- function(adjlist, directed, selfloops, weighted=NULL){
  if(directed)
    mode <- 'directed'
  if(!directed)
    mode <- 'undirected'

  lapply(X = adjlist, FUN = igraph::graph_from_adjacency_matrix, mode=mode, diag=selfloops, weighted=weighted)
}


#' @describeIn ghype Fitting ghype models from an igraph graph
#'
#' @export
#'
#'
ghype.igraph <- function(graph, directed, selfloops, xi=NULL, omega=NULL, unbiased=FALSE, regular=FALSE, ...){
  if(igraph::is_bipartite(graph)){
    adj <- igraph::get.incidence(graph = graph, sparse = FALSE)
  } else{
    adj <- igraph::get.adjacency(graph = graph, type = "upper", sparse = FALSE)
    if(!directed)
      adj <- adj + t(adj)
  }

  if(is.null(xi)){
    xi=ComputeXi(adj,directed,selfloops)
  }

  if(nrow(adj)==ncol(adj)){
    n <- nrow(adj)
  } else{
    n <- nrow(adj)+ncol(adj)
  }
  df <- regular + (!regular) * (1+directed) * n

  if(is.null(omega)){
    if(unbiased){
      omega <- matrix(1,nrow(adj), ncol(adj))
    } else{
      omega <- FitOmega(adj = adj, xi = xi, directed = directed, selfloops = selfloops)
      df <- df + sum(mat2vec.ix(omega,directed,selfloops))
    }
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
                         'regular' = regular,
                         'unbiased' = unbiased,
                         'df' = df), ...)
  return(model)
}


#' BootstrapProperty computes igraph analytics function on ensemble
#'
#' @param graph igraph graph
#' @param property igraph function that can be applied to a graph
#' @param directed boolean
#' @param selfloops boolean
#' @param nsamples number of samples from ensemble. defaults to 1000
#' @param xi matrix, default null
#' @param omega matrix, default null
#' @param model ghype model from which to extract xi and omega, default to null
#' @param m int, number of edges to sample from model
#' @param seed seed
#' @param ... other parameters to pass to `property`
#'
#' @return
#'
#' vector of length nsamples
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' library(igraph)
#' data('adj_karate')
#' result <- BootstrapProperty(adj_karate, page_rank, FALSE, FALSE, nsamples=100)
#' }
#'
BootstrapProperty <- function(graph, property, directed, 
            selfloops, nsamples=1000, xi=NULL, omega=NULL, 
            model=NULL, m=NULL, seed=NULL, ...){

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
    model <- ghype(graph = graph, directed, selfloops, xi, omega)

  rsamples <- rghype(nsamples, model, m, seed=seed)
  gsamples <- CreateIgGraphs(adjlist = rsamples, directed = directed, selfloops = selfloops)
  if(as.character(substitute(property)) %in% functionslist){
    dproperty <- sapply(X = gsamples, FUN = function(graph, directed, ...){match.fun(FUN = property)(graph, directed=directed, ...)$vector}, directed=directed, ...)
  } else{
    dproperty <- sapply(X = gsamples, FUN = property, ...)
  }
  return(dproperty)
}
