#' Fit propensity matrix for full model
#'
#' (auxiliary function)
#'
#' @param adj adjacency matrix
#' @param xi combinatorial matrix
#' @param directed boolean
#' @param selfloops boolean
#'
#' @return
#' propensity matrix
#'
#' @export
#'
#' @examples 
#' data(adj_karate)
#' xi <- compute_xi(adj_karate, FALSE, FALSE)
#' FitOmega(adj_karate, xi, FALSE, FALSE)
#' 
FitOmega <- function(adj, xi, directed, selfloops){
  return(MLE_omega_wallenius(adj,xi,directed,selfloops))
}


MLE_omega_idx <- function(adj,xi){
  ## provides vector of indices of elements of adj for which omega needs computing.
  ## omega for adj==Xi has to be set to 1 (no perfect fit)
  ## omega for adj==0 has to be set to 0 to speed up process
  ## if Xi==0 and adj!=0 rise warning and set omega to 0 anyway (MLE does not exist)
  idx.zero <- adj==0 | xi==0
  idx.one <- adj>=xi
  if(any(adj[xi<=adj]!=0))
    warning('MLE does not exist -- adj[xi<=adj]!=0')
  return(list("zero"=idx.zero,"one"=idx.one))
}

MLE_omega_wallenius <- function(adj,xi,directed,selfloops){
  omega.matrix <- matrix(0,nrow(adj),ncol(adj))
  idx <- mat2vec.ix(adj,directed,selfloops)
  adj.v <- as.vector(adj[idx])
  xi.v <- as.vector(xi[idx])

  omega.v <- rep(NA,length(adj.v))
  tmp <- MLE_omega_idx(adj.v,xi.v)
  idx.zero <- tmp$zero; idx.one <- tmp$one; rm(tmp)
  omega.v[idx.one] <- 1; omega.v[idx.zero] <- 0
  omega.v[!idx.one & !idx.zero] <- fitted_omega_wallenius(adj.v[!idx.one & !idx.zero], xi.v[!idx.one & !idx.zero])

  omega.matrix[idx] <- omega.v

  return(omega.matrix)
}


fitted_omega_wallenius <- function(adj,xi){
  ### fits omega VECTOR given adj and xi VECTORS already cleaned
  a <- log(1 - adj / xi)
  k <- min(a, na.rm = T)
  omega <- a / k
}
