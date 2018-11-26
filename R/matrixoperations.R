#' Auxilliary function, gives mask for matrix for directed,
#' undirected etc.
#'
#' @param mat  matrix
#' @param directed  a boolean argument specifying whether object is directed or not.
#' @param selfloops  a boolean argument specifying whether the model should incorporate selfloops.
#' @return
#' a boolean matrix that can be used to mask adjacency matrices.
#'
#' @export
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


#' Auxilliary function, produces matrix from vector
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
vec2mat <- function(vec,directed,selfloops,n){
  if(length(n)>1){
    mat <- matrix(0,n[2],n[3])
  } else{
    mat <- matrix(0,n,n)
  }

  idx <- mat2vec.ix(mat,directed,selfloops)
  mat[idx] <- vec
  return(mat)
}
