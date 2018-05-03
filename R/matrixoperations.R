#' Auxilliary function, gives indices of adjacency matrix for directed,
#' undirected etc.
#'
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param mat  ~~Describe \code{mat} here~~
#' @param directed  ~~Describe \code{directed} here~~
#' @param selfloops  ~~Describe \code{selfloops} here~~
#' @return
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
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
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param vec  ~~Describe \code{mat} here~~
#' @param directed  ~~Describe \code{directed} here~~
#' @param selfloops  ~~Describe \code{selfloops} here~~
#' @return
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#'
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
