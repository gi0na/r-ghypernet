#' Compute Phi Matrix
#'
#' @param graph n x n adjacency matrix of the graph of which computing the phi matrix
#' @param model ghype model against which computing the phi matrix
#' @param lightMemory boolean, store model and graph? default FALSE
#'
#' @return n x n phi matrix
#' @export
#'
#' @examples
#' 
#' data("adj_karate")
#' model <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' (phi <- phi_matrix(graph = adj_karate, model = model))
#' plot(phi)
#' 
phi_matrix <- function(graph, model, lightMemory=FALSE){
  # compute under and over represented dyads
  over <- linkSignificance(graph = graph, model=model, under=FALSE)
  under <- linkSignificance(graph = graph, model=model, under=TRUE)
  
  #compute phi
  phi_mat <- under - over
  
  if(lightMemory){
    phi <- list(matrix=phi_mat, graph=NULL, model=NULL)
  } else{
    phi <- list(matrix=phi_mat, graph=graph, model=model)
  }
  class(phi) <- append('phi_matrix',class(phi))
  
  return(phi)
}


#' Print method for objects of class phi_matrix
#'
#' @param x an object of class phi_matrix
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
#' @examples
#' data("adj_karate")
#' model <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' phi <- phi_matrix(graph = adj_karate, model = model)
#' print(phi)
#' 
print.phi_matrix <- function(x, ...){
  rownames(phi$matrix) <- colnames(phi$matrix) <- rownames(phi$graph)
  print(phi$matrix)
}

#' Plot method for objects of class phi_matrix
#'
#' @param x an object of class phi_matrix
#' @param ... arguments to be passed to methods, such as graphical parameters
#'
#' @export
#' @example 
#' data("adj_karate")
#' model <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' phi <- phi_matrix(graph = adj_karate, model = model)
#' plot(phi)
#' 
plot.phi_matrix <- function(x, ...){
  rownames(phi$matrix) <- colnames(phi$matrix) <- rownames(phi$graph)
  hmcol<-RColorBrewer::brewer.pal(11,"RdBu")
  heatmap(phi$matrix, symm = TRUE, col=hmcol, ...)
}
