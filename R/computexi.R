#' Auxilliary function. Computes combinatorial matrix.
#'
#' Combinatorial matrix computed according to soft
#' configuration model or 'regular' gnp model.
#'
#' @param adj adjacency matrix
#' @param directed boolean
#' @param selfloops boolean
#' @param regular boolean. Is the combinatorial matrix computed for configuration model or for regular gnp model? default FALSE.
#'
#' @return
#' combinatorial matrix
#'
#' @export
#'
ComputeXi <- function(adj, directed, selfloops, regular = FALSE) {
  # returns the matrix xi according to the nodes degrees

  if(!directed & !isSymmetric(adj)){
    warning('Trying to compute undirected ensemble for asymmetric adjacency matrix.
              Adjacency matrix symmetrised as adj <- adj + t(adj)')
    adj <- adj + t(adj)
  }

  if(regular){
    Kin <- apply(adj, 2, sum)
    ix <- mat2vec.ix(adj, directed, selfloops)
    m <- sum(adj[ix])
    M <- m^2
    if(!directed) M <- 4*M
    xiregular <- matrix(M/sum(ix), nrow(adj), ncol(adj))
    xiregular[!ix] <- 0
    if(!selfloops) diag(xiregular) <- 0
    xi <- ceiling(xiregular)
  } else{
    if(!selfloops) diag(adj) <- 0
    Kin <- apply(adj, 2, sum)
    Kout <- apply(adj, 1, sum)
    xi <- tcrossprod(Kout, Kin)
  }
  if(nrow(adj)==ncol(adj)){
    if(!selfloops & directed){
      diagxi <- diag(xi)
      vbas <- floor( diagxi /(ncol(adj)-1))
      diagxir <- diagxi -vbas*(ncol(adj)-1)
      diag(xi) <- 0
      xi <- #floor(
        (xi+t(vxi(1:nrow(adj), diagxir, vbas, xi, adj))) #/(as.integer(!directed)+1))
      # if(!directed){
      # xi[upper.tri(xi,FALSE)] <- ceiling( apply( cbind(xi[upper.tri(xi,FALSE)], t( xi )[upper.tri(xi,FALSE)] ), 1, mean) )
      # xi[lower.tri(xi,FALSE)] <- t( xi )[lower.tri(xi,FALSE)]
      # }
    } else {
      if(!directed){
        xi <- xi + t(xi) - diag(diag(xi))
        if(!selfloops){
          ix <- mat2vec.ix(adj, directed, selfloops)
          sdiag <- sum(diag(xi))
          toadd <- ceiling(sdiag/sum(Kin)*Kin/(nrow(adj)-1))
          # TODO: to improve. ugly for loop
          for(i in 1:nrow(xi)){
            xi[i,] <- xi[,i] <- xi[i,] + toadd[i]
          }
          diag(xi) <- 0
        }
      }
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
