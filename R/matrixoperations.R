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
