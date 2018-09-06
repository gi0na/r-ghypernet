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
fitBlockModel <- function(adj, labels, directed, selfloops, directedBlocks = FALSE, homophily = FALSE, inBlockOnly = FALSE, xi = NULL){

  if(!directed & (homophily | inBlockOnly) & directedBlocks)
    directedBlocks <- FALSE

  # generate unique blockids
  blockids <- as.numeric(
    plyr::mapvalues(labels, from = levels(factor(labels)),
                    to = c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))]))

  # generate block matrix
  blocks <- blockids %*% t(blockids)

  if(homophily & inBlockOnly){
    e <- simpleError('homophily and inBlockOnly parameters are in conflict')
    stop(e)
  }

  # if homophily model, fit only in-group vs out-group
  if(homophily){
    blocks[blocks %in% blockids^2] <- 1
    blocks[blocks != 1] <- 2
  }

  # if inBlockOnly model, fit different in-group parameters
  # but only one out-group parameter
  if(inBlockOnly){
    blocks[!(blocks %in% blockids^2)] <- 0
  }

  if(directedBlocks){
    blocks[lower.tri(blocks,F)] <- - blocks[lower.tri(blocks,F)]
    blocks[abs(blocks) %in% unique(blockids)^2] <- abs(blocks[abs(blocks) %in% unique(blockids)^2])
  }

  # construct xi matrix
  if(is.null(xi)){
    xi <- ComputeXi(adj,directed,selfloops)
  }

  # generate map xi and adj values to blocks
  xiframe <- data.frame(xi=xi[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])
  adjframe <- data.frame(adj=adj[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])

  # compute xi-blocks and adj-blocks
  xib <- plyr::ddply(xiframe, 'block', plyr::numcolwise(sum))
  mb <- plyr::ddply(adjframe, 'block', plyr::numcolwise(sum))
  m <- sum(adj[mat2vec.ix(xi,directed,selfloops)])

  # solve linear system for omega:
  # build matrix and vector

  # A <- - mb[,2] %*% t(xib[,2]) + diag(m*xib[,2])
  #
  # # remove entries = 0 to set omegas=0 there
  # zerosid <- which(mb[,2]==0)
  #
  # if(length(zerosid)!=0){
  #   A <- A[-zerosid,]
  #   A <- A[,-zerosid]
  # }
  #
  #
  # # set one parameter omega to 1
  # b <- - A[-1,1]
  #
  #
  # # solve system to find omega values
  # omegab <- rep(0,length(xib[,1]))
  # if(length(zerosid)!=0){
  #   omegab[-zerosid] <- c(1,solve(A[-1,-1],b))
  # } else{
  #   omegab <- c(1,solve(A[-1,-1],b))
  # }

  # use ensemble for omega:

  omega.v <- rep(NA,length(mb[,2]))
  tmp <- MLE_omega_idx(mb[,2],xib[,2])
  idx.zero <- tmp$zero; idx.one <- tmp$one; rm(tmp)
  omega.v[idx.one] <- 1; omega.v[idx.zero] <- 0
  omega.v[!idx.one & !idx.zero] <- fitted_omega_wallenius(mb[,2][!idx.one & !idx.zero], xib[,2][!idx.one & !idx.zero])

  # omegab <- FitOmega(adj = vec2mat(vec = mb[,2],directed = directedBlocks,selfloops = TRUE,(length(unique(blockids))*(!homophily)+2*homophily)),
  #          xi = vec2mat(vec = xib[,2],directed = directedBlocks,selfloops = TRUE,(length(unique(blockids))*(!homophily)+2*homophily)),
  #          directed = directedBlocks, selfloops = TRUE)
  # omegab <- omegab[mat2vec.ix(mat = omegab, directed = directedBlocks, selfloops = TRUE)]
  #
  # if(homophily)
  #   omegab <- omegab[-3]

  omegab <- data.frame(block=xib[,1],omegab=omega.v)

  # map values to full omega vector
  omegav <- plyr::mapvalues(xiframe[,2],from=sort(unique(xiframe[,2])), to=omegab[,2][rank(omegab[,1])])

  # generate omega matrix
  omega <- vec2mat(omegav,directed,selfloops,nrow(adj))

  # generate and return ensemble
  model <- ghype(object = adj, directed = directed, selfloops = selfloops, xi = xi, omega = omega)

  # generate block omega matrix for reference
  if((!homophily) & (!inBlockOnly)){
    blockOmega <- c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))] %*% t(c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))])

    if(directedBlocks)
      blockOmega[lower.tri(blockOmega,F)] <- - blockOmega[lower.tri(blockOmega,F)]

    rownames(blockOmega) <- colnames(blockOmega) <- levels(factor(labels))
    blockOmega <- plyr::mapvalues(blockOmega,from=unique(sort(blockOmega)), to=omegab[,2][rank(omegab[,1])])
  } else{
    blockOmega <- NULL
  }
  model$blockOmega <- blockOmega
  model$df <- sum(mat2vec.ix(xi,directed,selfloops)) + nrow(omegab)-1
  model$directedBlocks <- directedBlocks
  ci <- cbind(rep(0,length(xib[,1])),rep(0,length(xib[,1])),rep(0,length(xib[,1])))
  # if(length(zerosid)!=0){
  #   ci[-zerosid,][1,] <- c(1,1,0)
  #   ci[-zerosid,][-1,] <-
  #     blockmodel.ci(omegaBlocks = omegab[-zerosid,2][-1],
  #                   xiBlocks = xib[-zerosid,2][-1],
  #                   mBlocks = mb[-zerosid,2][-1], m=model$m)
  # } else{
    ci[1,] <- c(1,1,0)
    ci[-1,] <-
      blockmodel.ci(omegaBlocks = omegab[,2][-1],
                    xiBlocks = xib[,2][-1],
                    mBlocks = mb[,2][-1], m=model$m)
  # }
  model$ci <- ci
  model$coef <- omegab[,2]

  class(model) <- append('ghypeBlock', class(model))
  model$call <- match.call()

  return(model)
}




#' Computes Fisher Information matrix for estimators in block models.
#'
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param beta  ~~Describe \code{beta} here~~
#' @param w  ~~Describe \code{w} here~~
#' @param xi  ~~Describe \code{xi} here~~
#' @param adj  ~~Describe \code{adj} here~~
#' @param directed  ~~Describe \code{directed} here~~
#' @param selfloops  ~~Describe \code{selfloops} here~~
#' @return val
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#'
#'
JnBlock <- function(omegaBlocks, xiBlocks, mBlocks, m) {
  # Returns Fisher Information
  if(length(xiBlocks)>1){
    Jn <- (m*diag(xiBlocks) - mBlocks %*% t(xiBlocks))
  } else{
    Jn <- (m*xiBlocks - mBlocks %*% t(xiBlocks))
  }
  return(Jn)
}


#' Confidence intervals for block models.
#'
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param nr.m  ~~Describe \code{nr.m} here~~
#' @param w  ~~Describe \code{w} here~~
#' @param adj  ~~Describe \code{adj} here~~
#' @param pval  ~~Describe \code{pval} here~~
#' @return
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#'
blockmodel.ci <- function(omegaBlocks, xiBlocks, mBlocks, m,
                  pval=.05) {
  jn <- JnBlock(omegaBlocks, xiBlocks, mBlocks, m)
  jn <- sqrt(diag(solve(jn)))

  ci <- cbind(omegaBlocks - stats::qnorm(pval/2,
                                  lower.tail = F) * jn, omegaBlocks +
                stats::qnorm(pval/2, lower.tail = F) *
                jn, jn)
  colnames(ci) <- c(paste(pval, "% ci"), paste(1 - pval, "% ci"), "jn")
  ci
}
