#' Compute Phi Elementwise with hyperg distr
#'
#' @param edgecount 
#' @param xi 
#' @param sum_xi 
#' @param m 
#'
#' @return
#'
phi_element_hyperg <- function(edgecount, xi, sum_xi, m){
  2 * phyper(q = edgecount, m = xi, n = sum_xi-xi, k = m) - dhyper(x = edgecount, m = xi, n = sum_xi-xi, k = m) - 1
}

#' Compute Phi Elementwise with binomial distr
#'
#' @param edgecount 
#' @param p 
#' @param m 
#'
#' @return
#'
phi_element_binomial <- function(edgecount, p, m){
  2 * pbinom(q = edgecount, size = m, prob = p) - dbinom(x = edgecount, size = m, prob = p) - 1
}

#' Compute Phi Elementwise with wallenius distr
#'
#' @param edgecount 
#' @param xi 
#' @param sum_xi 
#' @param oddsratio 
#' @param m 
#'
#' @return
#'
phi_element_wallenius <- function(edgecount, xi, sum_xi, oddsratio, m){
  2 * BiasedUrn::pWNCHypergeo(x = edgecount, m1 = xi, m2 = sum_xi-xi, n = m, odds = oddsratio) - 
    BiasedUrn::dWNCHypergeo(x = edgecount, m1 = xi, m2 = sum_xi-xi, n = m, odds = oddsratio) - 1
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
#' data("adj_karate")
#' model <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' (phi <- phi_matrix(graph = adj_karate, model = model))
#' hmcol<-RColorBrewer::brewer.pal(11,"RdBu")
#' plot(phi, col=hmcol, legend="col")
phi_matrix <- function(graph, model, lightMemory=FALSE, useBinomial=FALSE){
  
  idx <- mat2vec.ix(graph, model$directed, model$selfloops)
  
  # compute parameters for marginal distributions
  sum_xi <- sum(model$xi[idx])
  m <- model$m
  
  # choose distribution
  useHyperg <- FALSE
  useWallenius <- FALSE
  if(!useBinomial){
    if(all(model$omega[idx]==model$omega[idx][1])){
      useHyperg <- TRUE
    } else{
      if( requireNamespace("BiasedUrn", quietly = TRUE) && (sum_xi-mean(model$xi[idx]))/m<1e4 )
        useWallenius <- TRUE
    }
  }
  if(!useHyperg & !useWallenius)
    useBinomial <- TRUE
  
  # compute values using the correct distr
  if(useHyperg){
    # compute under and over represented dyads
    phi_vals <- Vectorize(phi_element_hyperg, vectorize.args = c('edgecount', 'xi'))(
      edgecount = graph[idx], xi = model$xi[idx], sum_xi = sum_xi, m = model$m)
  } else{
    #compute oddsratio vector
    xibar <- sum_xi-model$xi[idx]
    omegabar <- (sum(model$xi[idx]*model$omega[idx])-model$xi[idx]*model$omega[idx])/xibar
    
    if(useWallenius){
      # compute under and over represented dyads
      phi_vals <- Vectorize(phi_element_wallenius, vectorize.args = c('edgecount', 'xi', 'oddsratio'))(
        edgecount = graph[idx], xi = model$xi[idx], sum_xi = sum_xi, oddsratio = model$omega[idx]/omegabar, m = model$m)
    
    } else{ #use binomial
      # get p vector
      p <- model$xi[idx]*model$omega[idx] / (
        model$xi[idx] * model$omega[idx]+xibar*omegabar
      )
      # compute under and over represented dyads
      phi_vals <- Vectorize(phi_element_binomial, vectorize.args = c('edgecount', 'p'))(
        edgecount = graph[idx], p = p, m = model$m)
    }
  }
  
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
  rownames(x$matrix) <- colnames(x$matrix) <- rownames(x$graph)
  print(x$matrix)
}

#' Plot method for objects of class phi_matrix
#'
#' @param x an object of class phi_matrix
#' @param ... arguments to be passed to methods, such as graphical parameters
#'
#' @export
#' @examples
#' data("adj_karate")
#' model <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' phi <- phi_matrix(graph = adj_karate, model = model)
#' hmcol<-RColorBrewer::brewer.pal(11,"RdBu")
#' plot(phi, col=hmcol)
#' 
plot.phi_matrix <- function(x, ...){
  rownames(phi$matrix) <- colnames(phi$matrix) <- rownames(phi$graph)
  heatmap(phi$matrix, symm = TRUE, ...)
}
