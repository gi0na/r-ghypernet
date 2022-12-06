#' Fitting bccm models
#'
#' bccm is used to fit a block-constrained configuration model.
#'
#'
#' @param adj the adjacency matrix of the graph.
#' @param labels vector or list. contains the vertex labels to generate the blocks in the bccm. In the case of bipartite graphs should be a list of two vectors, the first one with row labels and the second one with column labels.
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
  # check if labels are all identicals
  specs <- check_specs.matrix(adj)
  bipartite <- specs[3]
  if(is.null(directed) | is.null(selfloops)){
    if(is.null(directed)) directed <- specs[1]
    if(is.null(selfloops)) selfloops <- specs[2]
  }
  if(bipartite){
    directed <- specs[1]
    selfloops <- specs[2]
  }
  
  if(length(unique(as.vector(labels))) == 1){
    # return unbiased ghype
    model <- ghype.matrix(graph = adj, directed = directed, selfloops = selfloops, xi = xi, unbiased = TRUE, regular = regular, ...)
    return(model)
  }
  model <- .bccm(adj = adj, labels = labels, directed = directed, selfloops = selfloops, directedBlocks = directedBlocks, homophily = homophily, inBlockOnly = inBlockOnly, xi = xi, regular = regular, bipartite = bipartite, ...)
  model$call <- match.call()
  return(model)
}

.bccm <- function(adj, labels, directed, selfloops, directedBlocks, homophily, inBlockOnly, xi, regular, ignore_pvals=FALSE, bipartite, ...){

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
  if(!bipartite){
    blockids <- as.numeric(
      plyr::mapvalues(labels, from = levels(factor(labels)),
                      to = c(1,numbers::Primes(length(labels)*8))[1:(length(unique(labels)))]))
  
    # generate block matrix
    blocks <- blockids %*% t(blockids)
    
    labs_map <- tibble(labs=labels,
                       ids=blockids, matchme=1
    ) %>%
      group_by(.data$labs, .data$ids) %>% 
      summarise(matchme = 1, .groups = 'drop')
    
  } else{
    blockids1 <- as.numeric(plyr::mapvalues(labels[[1]], from = levels(factor(unlist(labels))), 
                                            to = c(1, numbers::Primes(length(unlist(labels)) * 8))[1:(length(unique(unlist(labels))))], warn_missing = FALSE))
    blockids2 <- as.numeric(plyr::mapvalues(labels[[2]], from = levels(factor(unlist(labels))), 
                                            to = c(1, numbers::Primes(length(unlist(labels)) * 8))[1:(length(unique(unlist(labels))))], warn_missing = FALSE))
    
    blocks <- blockids1 %*% t(blockids2)
    
    labs_map1 <- tibble(labs=labels[[1]],
                       ids=blockids1, matchme=1
    ) %>%
      group_by(.data$labs, .data$ids) %>% 
      summarise(matchme = 1, .groups = 'drop')
    
    labs_map2 <- tibble(labs=labels[[2]],
                        ids=blockids2, matchme=1
    ) %>%
      group_by(.data$labs, .data$ids) %>% 
      summarise(matchme = 1, .groups = 'drop')
  }
  
  if(directedBlocks){
    full_join(labs_map,labs_map,  by = 'matchme') %>% select(-"matchme") %>%
      rowwise() %>%
      mutate(lab = paste(.data$labs.x,.data$labs.y, sep = '->'),
             id = if_else(.data$ids.x<=.data$ids.y,.data$ids.x*.data$ids.y,-.data$ids.x*.data$ids.y)) %>%
      select("labs.x","labs.y","lab","id") %>%
      ungroup() ->
      labs_map

    blocks[order(blockids),order(blockids)][lower.tri(blocks,F)] <- -blocks[order(blockids),order(blockids)][lower.tri(blocks,F)]
    # blocks[lower.tri(blocks,F)] <- - blocks[lower.tri(blocks,F)]
    blocks[abs(blocks) %in% unique(blockids)^2] <- abs(blocks[abs(blocks) %in% unique(blockids)^2])
  } else{
    if(!bipartite){
      full_join(labs_map,labs_map,  by = 'matchme') %>% select(-"matchme") %>%
        rowwise() %>%
        mutate(lab = paste(.data$labs.x,.data$labs.y, sep = '<->'),
               id = .data$ids.x*.data$ids.y) %>%
        group_by(id) %>%
        summarise(labs.x = first(.data$labs.x),labs.y = first(.data$labs.y), lab = first(.data$lab)) ->
        labs_map
    }
    if(bipartite){
      full_join(labs_map1,labs_map2,  by = 'matchme') %>% select(-"matchme") %>%
        rowwise() %>%
        mutate(lab = paste(.data$labs.x,.data$labs.y, sep = '<->'),
               id = .data$ids.x*.data$ids.y) %>%
        group_by(id) %>%
        summarise(labs.x = first(.data$labs.x),labs.y = first(.data$labs.y), lab = first(.data$lab)) ->
        labs_map
    }
  }

  if(homophily & inBlockOnly){
    e <- simpleError('homophily and inBlockOnly parameters are in conflict')
    stop(e)
  }

  # if homophily model, fit only in-group vs out-group
  if(homophily){
    if(!bipartite){
      blocks[blocks %in% unique(blockids)^2] <- 1
      blocks[blocks != 1] <- 2
      labs_map %>%
        mutate(homo = id %in% unique(blockids)^2) %>%
        mutate(id = (!.data$homo) + 1) %>%
        group_by(id) %>%
        summarise(labs.x = first(.data$labs.x),labs.y = first(.data$labs.y), lab = if_else(first(id) == 1,'homologue', 'hetero')) ->
        labs_map
    } else{
      blocks[blocks %in% unique(c(blockids1,blockids2))^2] <- 1
      blocks[blocks != 1] <- 2
      labs_map %>%
        mutate(homo = id %in% unique(c(blockids1,blockids2))^2) %>%
        mutate(id = (!.data$homo) + 1) %>%
        group_by(id) %>%
        summarise(labs.x = first(.data$labs.x),labs.y = first(.data$labs.y), lab = if_else(first(id) == 1,'homologue', 'hetero')) ->
        labs_map
    }
  }

  # if inBlockOnly model, fit different in-group parameters
  # but only one out-group parameter
  if(inBlockOnly){
    if(!bipartite){
      blocks[!(blocks %in% blockids^2)] <- 0
      labs_map %>%
        mutate(inblock = id %in% unique(blockids)^2,
               id = if_else(.data$inblock,id,0),
               lab = if_else(.data$inblock,.data$lab,'betweenblocks')) %>%
        group_by(id) %>%
        summarise(labs.x = first(.data$labs.x),labs.y = first(.data$labs.y), lab = first(.data$lab),id = first(id)) ->
        labs_map
    } else{
      blocks[!(blocks %in% unique(c(blockids1,blockids2))^2)] <- 0
      labs_map %>%
        mutate(inblock = id %in% unique(c(blockids1,blockids2))^2,
               id = if_else(.data$inblock,id,0),
               lab = if_else(.data$inblock,.data$lab,'betweenblocks')) %>%
        group_by(id) %>%
        summarise(labs.x = first(.data$labs.x),labs.y = first(.data$labs.y), lab = first(.data$lab),id = first(id)) ->
        labs_map
    }
  }
  # construct xi matrix
  if(is.null(xi)){
    if(is.null(regular)){
      regular <- FALSE
      if(conf.test(graph = adj, directed = directed, selfloops = selfloops)$p.value>1e-4)
        regular <- TRUE
    }

    xi <- compute_xi(adj,directed,selfloops, regular=regular)

  } else{
    if(length(xi) == 1 && xi == 'regular'){
      k1 = rep(sum(adj)/nrow(adj),nrow(adj))
      k2 = rep(sum(adj)/ncol(adj),ncol(adj))
      xi <- round(k1 %*% t(k2))
    }
  }

  # generate map xi and adj values to blocks
  xiframe <- data.frame(xi=xi[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])
  adjframe <- data.frame(adj=adj[mat2vec.ix(xi,directed,selfloops)], block=blocks[mat2vec.ix(xi,directed,selfloops)])

  # compute xi-blocks and adj-blocks
  xib <- xiframe %>% group_by(.data$block) %>% summarise(xi=sum(xi))
  mb <- adjframe %>% group_by(.data$block) %>% summarise(adj=sum(adj))

  ## BUG: if one node singleton in community and selfloops not allowed there is no entry in omegab for it, set manually to 0
  if(!selfloops & !bipartite){
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


  omegab <- data.frame(block=xib$block,omegab=omega.v) %>%
    left_join(labs_map, by = c('block'='id'))

  # map values to full omega vector
  omegav <- left_join(xiframe, omegab %>% select('block','omegab'), by = 'block') %>% pull(omegab)
  # omegav <- plyr::mapvalues(xiframe$block,from=sort(xiframe$block), to=omegab$omega[order(omegab$block)])

  # generate omega matrix
  if(bipartite){ # if bipartite
    n <- c(nrow(adj)+ncol(adj),nrow(adj),ncol(adj))
  } else{
    n <- nrow(adj)
  }
  omega <- vec2mat(omegav,directed,selfloops,n)

  # generate and return ensemble
  model <- ghype(graph = adj, directed = directed, selfloops = selfloops, xi = xi, omega = omega, regular = regular, ...)

  # generate block omega matrix for reference
  if( (!homophily) & (!inBlockOnly)){
    if(bipartite){ # if bipartite
      blockOmega <- sort(unique(blockids1)) %*% t(sort(unique(blockids2)))
      rownames(blockOmega) <- levels(factor(labels[[1]]))
      colnames(blockOmega) <- levels(factor(labels[[2]]))
      blockOmega <- plyr::mapvalues(blockOmega, from = unique(sort(blockOmega)), 
                                    to = omegab[, 2][rank(omegab[, 1])])
    } else{
      
      omegab %>% arrange(.data$labs.y,.data$labs.x) %>% 
        pull(omegab) %>%
        vec2mat(directed = directedBlocks, selfloops = TRUE, n = length(unique(labels))) ->
        blockOmega
      if(!directedBlocks){
        diag(blockOmega) <- diag(blockOmega)/2
      }
      
      rownames(blockOmega) <- colnames(blockOmega) <- sort(unique(labels))
    }
  } else{
    blockOmega <- NULL
  }
  model$blockOmega <- blockOmega
  # model$df <- regular + (1-regular)*nrow(xi)*(1+directed) + nrow(omegab)-1 - sum(xib[,2]==0)
  model$df <- regular + (1 - regular) * (nrow(xi) + directed*ncol(xi)) + nrow(omegab) - 1 - sum(xib[, 2] == 0)
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
    ci[1,] <- c(omegab[1,'omegab'],omegab[1,'omegab'],0)
    ci[-1,] <-
      blockmodel.ci(omegaBlocks = omegab$omegab[-1],
                    xiBlocks = xib$xi[-1],
                    mBlocks = mb$adj[-1], m=model$m)
  }
  # }
  model$ci <- ci
  model$coef <- omegab[,'omegab']
  names(model$coef) <- omegab[,'lab']
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
