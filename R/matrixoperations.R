#' Auxiliary function, gives mask for matrix for directed,
#' undirected etc.
#'
#' @param mat  matrix
#' @param directed  a boolean argument specifying whether object is directed or not.
#' @param selfloops  a boolean argument specifying whether the model should incorporate selfloops.
#' @return
#' a boolean matrix that can be used to mask adjacency matrices.
#'
#' @export
#' 
#' @examples
#' data('adj_karate')
#' mat2vec.ix(adj_karate, FALSE, FALSE)
#' 
mat2vec.ix <- function(mat, directed,
                       selfloops) {
  # Returns the indices to
  # vectorise adjacency matrices
  # removing unused entries in the
  # case of undirected or no
  # selfloops graphs
  if(nrow(mat)==ncol(mat)){
    if (!directed) {
      mat[lower.tri(mat, !selfloops)] <- NA
    } else {
      if (!selfloops)
        diag(mat) <- NA
    }
  }
  ix <- !is.na(mat)
  return(ix)
}


#' Auxiliary function, produces matrix from vector
#' 
#' The number of elements of vec are the number of non-zero elements in the
#' adjacency matrix.
#' It performs the opposite operation of `mat2vec.ix`.
#'
#'
#' @param vec  vector to be put in matrix form
#' @param directed  a boolean argument specifying whether object is directed or not.
#' @param selfloops  a boolean argument specifying whether the model should
#' incorporate selfloops.
#' @param n vector. if length(n)==1, n is the number of vertices. If length(n)==3
#' first element is number of vertices, second and third elements are number of
#' vertices for row and column of bipartite matrix.
#' @return
#' matrix nxn generated from vector.
#' @export
#' 
#' @examples
#' data('adj_karate')
#' ix <- mat2vec.ix(adj_karate, FALSE, FALSE)
#' vec <- adj_karate[ix]
#' vec2mat(vec, FALSE, FALSE, nrow(adj_karate))
#' 
vec2mat <- function(vec,directed,selfloops,n){
  if(length(n)>1){
    mat <- matrix(0,n[2],n[3])
  } else{
    mat <- matrix(0,n,n)
  }

  idx <- mat2vec.ix(mat,directed,selfloops)
  mat[idx] <- vec
  if(!directed) mat <- mat + t(mat)
  return(mat)
}

#' Maps adjacency matrix to edgelist
#'
#' @param adj matrix, the adjacency matrix
#' @param directed boolean, is the graph directed?
#'
#' @return a dataframe containing the edgelist
#' 
#' @examples 
#' data(contacts.adj)
#' el <- adj2el(contacts.adj)
#' 
#' @export
#' @import dplyr
#' @importFrom rlang .data
#'
adj2el <- function(adj, directed=TRUE){
  if(!directed) adj[lower.tri(adj,F)] <- 0
  reshape2::melt(adj, value.name = "value") %>%
    filter(.data$value!=0) %>%
    rename(sender='Var1', target='Var2', edgecount='value') %>%
    mutate_at(c("sender", "target"), as.character) -> el
  return(el)
}


#' Maps edgelist to adjacency matrix
#'
#' @param el dataframe containing a (weighted) edgelist. Column 1 is the sender, column 2 is the receiver, column 3 the number of edges.
#' @param nodes optional vector containing all node names in case disconnected nodes should be included.
#'
#' @return the (weighted) adjacency matrix corresponding the edgelist passed
#' @export
#'
el2adj <- function(el, nodes=NULL){
  if(is.null(nodes))
    nodes <- unique(c(el[,1], el[,2]))
  
  adj <- matrix(0, nrow = length(nodes), ncol = length(nodes),
                dimnames = list(nodes,nodes))
  adj[as.matrix(el[,1:2])] <- el[,3]
  return(adj)
}
