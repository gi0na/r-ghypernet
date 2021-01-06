#' Fitting bccm models
#'
#' bccm is used to fit a block-constrained configuration model.
#'
#'
#' @param adj the adjacency matrix of the graph.
#' @param labels vector, the vertex labels to generate the blocks in the bccm.
#' @param directed a boolean argument specifying whether the graph is directed or not.
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
#' @examples
#' data("vertexlabels","adj_karate")
#' blockmodel <- bccm(adj = adj_karate, labels = vertexlabels, directed = FALSE, selfloops = FALSE)
#' 
bccm <- function(adj, labels, directed = NULL, selfloops = NULL, directedBlocks = FALSE, homophily = FALSE, inBlockOnly = FALSE, xi = NULL, regular = FALSE, ...){
  model <- .bccm(adj = adj, labels = labels, directed = directed, selfloops = selfloops, directedBlocks = directedBlocks, homophily = homophily, inBlockOnly = inBlockOnly, xi = xi, regular = regular, ...)
  model$call <- match.call()
  return(model)
}

.bccm <- function(adj, labels, directed, selfloops, directedBlocks, homophily, inBlockOnly, xi, regular, ignore_pvals=FALSE, ...){
    if(is.null(directed) | is.null(selfloops)){
    specs <- check_specs(adj)
    if(is.null(directed)) directed <- specs[1]
    if(is.null(selfloops)) selfloops <- specs[2]
  }

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
  xib <- xiframe %>% group_by(.data$block) %>% summarise(xi=sum(xi))
  mb <- adjframe %>% group_by(.data$block) %>% summarise(adj=sum(adj))

  ## BUG: if one node singleton in community and selfloops not allowed there is no entry in omegab for it, set manually to 0
  if(!selfloops){
    vblockcounts <- plyr::count(blockids)
    if(any(vblockcounts$freq<2)){
      xib <- rbind(xib, cbind(block=vblockcounts$x[which(vblockcounts$freq<2)]^2, xi=rep(0, sum(vblockcounts$freq<2))))
      mb <- rbind(mb, cbind(block=vblockcounts$x[which(vblockcounts$freq<2)]^2, adj=rep(0, sum(vblockcounts$freq<2))))
    }
  }
  m <- sum(adj[mat2vec.ix(xi,directed,selfloops)])

  # use ensemble for omega:

  omega.v <- rep(NA,length(mb[,2]))
  tmp <- MLE_omega_idx(mb[,2],xib[,2])
  idx.zero <- tmp$zero; idx.one <- tmp$one; rm(tmp)
  omega.v[idx.one] <- 1; omega.v[idx.zero] <- 0

  omega.v[!idx.one & !idx.zero] <- fitted_omega_wallenius(mb[,2][!idx.one & !idx.zero], xib[,2][!idx.one & !idx.zero])



  omegab <- data.frame(block=xib$block,omegab=omega.v)

  # map values to full omega vector
  omegav <- plyr::mapvalues(xiframe$block,from=sort(unique(xiframe$block)), to=omegab$omega[rank(omegab$block[1:length(unique(xiframe$block))])])

  # generate omega matrix
  omega <- vec2mat(omegav,directed,selfloops,nrow(adj))

  # generate and return ensemble
  model <- ghype(graph = adj, directed = directed, selfloops = selfloops, xi = xi, omega = omega, regular = regular, ...)

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
  model$df <- regular + (1-regular)*nrow(xi)*(1+directed) + nrow(omegab)-1 - sum(xib[,2]==0)
  model$directedBlocks <- directedBlocks
  model$homophily <-  homophily
  model$inBlockOnly <- inBlockOnly
  ci <- cbind(rep(0,nrow(xib)),rep(0,nrow(xib)),rep(0,nrow(xib)))
  # if(length(zerosid)!=0){
  #   ci[-zerosid,][1,] <- c(1,1,0)
  #   ci[-zerosid,][-1,] <-
  #     blockmodel.ci(omegaBlocks = omegab[-zerosid,2][-1],
  #                   xiBlocks = xib[-zerosid,2][-1],
  #                   mBlocks = mb[-zerosid,2][-1], m=model$m)
  # } else{
  if(isFALSE(ignore_pvals)){
    ci[1,] <- c(omegab[1,2],omegab[1,2],0)
    ci[-1,] <-
      blockmodel.ci(omegaBlocks = omegab$omegab[-1],
                    xiBlocks = xib$xi[-1],
                    mBlocks = mb$adj[-1], m=model$m)
  }
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


# Confidence intervals for block models.
# Based on Fisher information matrix
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
