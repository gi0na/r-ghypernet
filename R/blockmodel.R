#' Fitting bccm models (deprecated)
#'
#' bccm is used to fit a block-constrained configuration model.
#'
#'
#' @param adj the adjacency matrix of the graph.
#' @param labels vector, the vertex labels to generate the blocks in the bccm.
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops boolean argument specifying whether the model should incorporate selfloops.
#' @param directedBlocks boolean argument specifying whether the model should incorporate directed blocks. Default to FALSE.
#' @param homophily boolean argument specifying whether the model should fit only homophily blocks. Default to FALSE.
#' @param inBlockOnly boolean argument specifying whether the model should fit only blocks over the diagonal. Default to FALSE.
#' @param xi an optional matrix defining the combinatorial matrix of the model.
#'
#' @return
#' bccm returns an object of class 'bccm' and 'ghype'.
#' 'bccm' objects expand 'ghype' objects incorporating the parameter estimates.
#' @export
#'
fitBlockModel <- function(adj, labels, directed, selfloops, directedBlocks = FALSE, homophily = FALSE, inBlockOnly = FALSE, xi = NULL){
  bccm(adj, labels, directed, selfloops, directedBlocks, homophily, inBlockOnly, xi)
}

#' Fitting bccm models
#'
#' bccm is used to fit a block-constrained configuration model.
#'
#'
#' @param adj the adjacency matrix of the graph.
#' @param labels vector, the vertex labels to generate the blocks in the bccm.
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops boolean argument specifying whether the model should incorporate selfloops.
#' @param directedBlocks boolean argument specifying whether the model should incorporate directed blocks. Default to FALSE.
#' @param homophily boolean argument specifying whether the model should fit only homophily blocks. Default to FALSE.
#' @param inBlockOnly boolean argument specifying whether the model should fit only blocks over the diagonal. Default to FALSE.
#' @param xi an optional matrix defining the combinatorial matrix of the model.
#' @param regular optional boolean, fit regular gnp model? if not specified chosen through lr.test.
#'
#' @return
#' bccm returns an object of class 'bccm' and 'ghype'.
#' 'bccm' objects expand 'ghype' objects incorporating the parameter estimates.
#' @export
#'
bccm <- function(adj, labels, directed = NULL, selfloops = NULL, directedBlocks = FALSE, homophily = FALSE, inBlockOnly = FALSE, xi = NULL, regular = NULL){

  specs <- check_specs(adj)
  directed <- specs[1]
  selfloops <- specs[2]

  if(is.matrix(adj)){
    if(!directed & !isSymmetric(adj)){
      warning('Trying to compute undirected ensemble for asymmetric adjacency matrix.
              Adjacency matrix symmetrised as adj <- adj + t(adj)')
      adj <- adj + t(adj)
    }
  }

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
    blocks[blocks %in% unique(blockids)^2] <- 1
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

  # # construct xi matrix
  # if(is.null(xi)){
  #   if(conf.test(graph = adj, directed = directed, selfloops = selfloops)>1e-4){
  #     ix <- mat2vec.ix(adj, directed, selfloops)
  #     m <- sum(adj[ix])
  #     xiregular <- matrix(m^2/sum(adj[ix]!=0), nrow(adj), ncol(adj))
  #     # if(!directed) xiregular <- xiregular + t(xiregular) - diag(diag(xiregular))
  #     xi <- ceiling(xiregular)
  #   } else{
  #     xi <- ComputeXi(adj,directed,selfloops)
  #   }
  # }

  # construct xi matrix
  if(is.null(xi)){
    if(is.null(regular)){
      regular <- FALSE
      if(conf.test(graph = adj, directed = directed, selfloops = selfloops)$p.value>1e-4)
        regular <- TRUE
    }

    xi <- ComputeXi(adj,directed,selfloops, regular=regular)

  } else{
    if(length(xi) == 1 && xi == 'regular')
      xi <- ComputeXi(adj,directed,selfloops, regular=TRUE)
  }

  # generate map xi and adj values to blocks
  xiframe <- data.frame(xi=xi[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])
  adjframe <- data.frame(adj=adj[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])

  # compute xi-blocks and adj-blocks
  xib <- plyr::ddply(xiframe, 'block', plyr::numcolwise(sum))
  mb <- plyr::ddply(adjframe, 'block', plyr::numcolwise(sum))

  ## BUG: if one node singleton in community and selfloops not allowed there is no entry in omegab for it, set manually to 0
  if(!selfloops){
    vblockcounts <- plyr::count(blockids)
    if(any(vblockcounts$freq<2)){
      xib <- rbind(xib, cbind(block=vblockcounts$x[which(vblockcounts$freq<2)]^2, xi=rep(0, sum(vblockcounts$freq<2))))
      mb <- rbind(mb, cbind(block=vblockcounts$x[which(vblockcounts$freq<2)]^2, adj=rep(0, sum(vblockcounts$freq<2))))
    }
  }
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
  # print(mb[,2][!idx.one & !idx.zero])
  # print(xib[,2][!idx.one & !idx.zero])
  # print(adjframe)
  # print(xiframe)
  # print('---')
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
  omegav <- plyr::mapvalues(xiframe[,2],from=sort(unique(xiframe[,2])), to=omegab[,2][rank(omegab[,1][1:length(unique(xiframe[,2]))])])

  # generate omega matrix
  omega <- vec2mat(omegav,directed,selfloops,nrow(adj))

  # generate and return ensemble
  model <- ghype(object = adj, directed = directed, selfloops = selfloops, xi = xi, omega = omega)

  # generate block omega matrix for reference
  if( (!homophily) & (!inBlockOnly)){
    blockOmega <- c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))] %*% t(c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))])

    if(directedBlocks)
      blockOmega[lower.tri(blockOmega,F)] <- - blockOmega[lower.tri(blockOmega,F)]

    rownames(blockOmega) <- colnames(blockOmega) <- levels(factor(labels))

    blockOmega <- plyr::mapvalues(blockOmega,from=unique(sort(blockOmega)), to=omegab[,2][rank(omegab[,1])])
  } else{
    blockOmega <- NULL
  }
  model$blockOmega <- blockOmega
  model$df <- xiregular + (1-xiregular)*nrow(xi)*(1+directed) + nrow(omegab)-1 - sum(xib[,2]==0)
  model$directedBlocks <- directedBlocks
  ci <- cbind(rep(0,length(xib[,1])),rep(0,length(xib[,1])),rep(0,length(xib[,1])))
  # if(length(zerosid)!=0){
  #   ci[-zerosid,][1,] <- c(1,1,0)
  #   ci[-zerosid,][-1,] <-
  #     blockmodel.ci(omegaBlocks = omegab[-zerosid,2][-1],
  #                   xiBlocks = xib[-zerosid,2][-1],
  #                   mBlocks = mb[-zerosid,2][-1], m=model$m)
  # } else{
    ci[1,] <- c(omegab[1,2],omegab[1,2],0)
    ci[-1,] <-
      blockmodel.ci(omegaBlocks = omegab[,2][-1],
                    xiBlocks = xib[,2][-1],
                    mBlocks = mb[,2][-1], m=model$m)
  # }
  model$ci <- ci
  model$coef <- omegab[,2]
  names(model$coef) <- omegab[,1]
  model$labels <- labels

  class(model) <- append(c('bccm','ghypeBlock'), class(model))
  model$call <- match.call()

  return(model)
}





#' Fisher Information matrix for estimators in block models.
#'
#' @param omegaBlocks the block parameters (vector)
#' @param xiBlocks the xi-block (vector)
#' @param mBlocks the adj-block (vector)
#' @param m the number of edges (scalar)
#'
#' @return
#' Fisher Information matrix
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
#' @return 3-vector
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#'
blockmodel.ci <- function(omegaBlocks, xiBlocks, mBlocks, m,
                  pval=.05) {
  jn <- JnBlock(omegaBlocks, xiBlocks, mBlocks, m)
  solve2 <- purrr::possibly(function(MM){
    sqrt(diag(solve(MM)))
  }, otherwise = NA)
  jn <- solve2(jn)

  ci <- cbind(omegaBlocks - stats::qnorm(pval/2,
                                  lower.tail = F) * jn, omegaBlocks +
                stats::qnorm(pval/2, lower.tail = F) *
                jn, jn)
  colnames(ci) <- c(paste(pval, "% ci"), paste(1 - pval, "% ci"), "jn")
  ci
}
