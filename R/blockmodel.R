#' Title
#'
#' @param adj
#' @param labels
#' @param directed
#' @param selfloops
#'
#' @return
#' @export
#'
#' @import numbers
#' @import plyr
#'
#' @examples
fitBlockModel <- function(adj, labels, directed, selfloops, directedBlocks = FALSE){

  if(!directed & directedBlocks)
    directedBlocks <- FALSE

  # generate unique blockids
  blockids <- as.numeric(
    plyr::mapvalues(labels, from = levels(factor(labels)),
                    to = c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))]))

  # generate block matrix
  blocks <- blockids %*% t(blockids)

  if(directedBlocks){
    blocks[lower.tri(blocks,F)] <- - blocks[lower.tri(blocks,F)]
    blocks[abs(blocks) %in% unique(blockids)^2] <- abs(blocks[abs(blocks) %in% unique(blockids)^2])
  }

  # construct xi matrix
  xi <- hypernets::ComputeXi(adj,directed,selfloops)

  # generate map xi and adj values to blocks
  xiframe <- data.frame(xi=xi[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])
  adjframe <- data.frame(adj=adj[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])

  # compute xi-blocks and adj-blocks
  xib <- plyr::ddply(xiframe, 'block', plyr::numcolwise(sum))
  mb <- plyr::ddply(adjframe, 'block', plyr::numcolwise(sum))
  m <- sum(adj[mat2vec.ix(xi,directed,selfloops)])

  # solve linear system for omega:
  # build matrix and vector
  A <- - mb[,2] %*% t(xib[,2]) + diag(m*xib[,2])
  b <- rep(0,nrow(A))

  # remove entries = 0 to set omegas=0 there
  zerosid <- which(mb[,2]==0)
  A <- A[-zerosid,]
  A <- A[,-zerosid]
  b <- b[-zerosid]

  # set max sized block to omega=1
  maxid <- which.max(mb[-zerosid,2])
  b[maxid] <- 1
  A[maxid,] <- b

  # solve system to find omega values
  omegab <- rep(0,length(xib[,1]))
  omegab[-zerosid] <- solve(A,b)
  omegab <- data.frame(block=xib[,1],omegab=omegab)

  # map values to full omega vector
  omegav <- plyr::mapvalues(xiframe[,2],from=sort(unique(xiframe[,2])), to=omegab[,2][rank(omegab[,1])])

  # generate omega matrix
  omega <- vec2mat(omegav,directed,selfloops,nrow(adj))

  # generate and return ensemble
  model <- ghype(object = adj, directed = directed, selfloops = selfloops, xi = xi, omega = omega)

  # generate block omega matrix for reference
  blockOmega <- c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))] %*% t(c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))])

  if(directedBlocks)
    blockOmega[lower.tri(blockOmega,F)] <- - blockOmega[lower.tri(blockOmega,F)]

  rownames(blockOmega) <- colnames(blockOmega) <- levels(factor(labels))
  blockOmega <- plyr::mapvalues(blockOmega,from=unique(sort(blockOmega)), to=omegab[,2][rank(omegab[,1])])
  model$blockOmega <- blockOmega
  model$df <- sum(mat2vec.ix(xi,directed,selfloops)) + nrow(omegab)-1

  return(model)
}
