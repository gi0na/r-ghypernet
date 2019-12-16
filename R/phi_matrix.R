#' Compute Phi Elementwise
#'
#' @param edgecount 
#' @param xi 
#' @param sum_xi 
#' @param m 
#'
#' @return
#' @export
#'
#' @examples
phi_element <- function(edgecount, xi, sum_xi, m){
  2 * phyper(q = edgecount, m = xi, n = sum_xi-xi, k = m) - dhyper(x = edgecount, m = xi, n = sum_xi-xi, k = m) - 1
}

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
#' 
phi_matrix <- function(graph, model, lightMemory=FALSE){
  # compute under and over represented dyads
  idx <- mat2vec.ix(graph, model$directed, model$selfloops)
  phi_vals <- Vectorize(phi_element, vectorize.args = c('edgecount', 'xi'))(
    edgecount = graph[idx], xi = model$xi[idx], sum_xi = sum(model$xi[idx]), m = model$m)
  phi_mat <- vec2mat(phi_vals, model$directed, model$selfloops, n = model$n)
  
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
#' hmcol<-RColorBrewer::brewer.pal(11,"RdBu")
#' plot(phi, col=hmcol, legend="col")
#' 
plot.phi_matrix <- function(x, ...){
  rownames(phi$matrix) <- colnames(phi$matrix) <- rownames(phi$graph)
  heatmap(phi$matrix, symm = TRUE, ...)
}
