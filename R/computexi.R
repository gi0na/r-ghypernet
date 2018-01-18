#' Auxilliary function. Computes Xi matrix.
#'
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param adj  ~~Describe \code{adj} here~~
#' @return
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' @export
ComputeXi <- function(adj, directed, selfloops) {
  # returns the matrix xi according to the nodes degrees
  Kin <- apply(adj, 2, sum)
  Kout <- apply(adj, 1, sum)
  xi <- tcrossprod(Kout, Kin)
  if(!selfloops & directed){
    diagxi <- diag(xi)
    vbas <- floor( diagxi /(ncol(adj)-1))
    diagxir <- diagxi -vbas*(ncol(adj)-1)
    set.seed(10)
    diag(xi) <- 0
    xi <- #floor(
      (xi+t(vxi(1:nrow(adj), diagxir, vbas, xi, adj))) #/(as.integer(!directed)+1))
    # if(!directed){
    # xi[upper.tri(xi,FALSE)] <- ceiling( apply( cbind(xi[upper.tri(xi,FALSE)], t( xi )[upper.tri(xi,FALSE)] ), 1, mean) )
    # xi[lower.tri(xi,FALSE)] <- t( xi )[lower.tri(xi,FALSE)]
    # }
  } else {
    if(!directed){
      xi <- ceiling(xi/2)
    }
  }
  return(xi)
}



vxi <- function(idx, diagxir, vbas, xi, adj){
  v <- rep(0, ncol(adj)-1)
  v[sample(ncol(adj)-1, diagxir[idx])] <- 1
  v <- v+vbas[idx]
  if(idx>1 & idx<ncol(xi)){
    v <- c(v[1:(idx-1)],0,v[idx:length(v)])
  } else{
    if(idx==1){
      v <- c(0,v)
    } else{
      v <- c(v,0)
    }
  }
  return(v)
}
vxi <- Vectorize(FUN = vxi, vectorize.args = 'idx')
