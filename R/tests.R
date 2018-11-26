#' Test regular (gnp) vs configuration model
#'
#' Likelihood ratio test for gnp vs configuration model.
#'
#' @param graph adjacency matrix or igraph graph
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param nempirical optional, number of graphs to sample from null distribution for empirical distribution.
#' @param parallel optional, number of cores to use or boolean for parallel computation.
#' If passed TRUE uses all cores-1, else uses the number of cores passed. If none passed
#' performed not in parallel.
#'
#' @return
#' p-value of test.
#'
#' @export
#'
conf.test <- function(graph, directed, selfloops, nempirical=NULL, parallel=NULL){

  adj <- graph
  if(requireNamespace("igraph", quietly = TRUE) && igraph::is.igraph(graph)){
    adj <- igraph::get.adjacency(graph, type='upper', sparse = FALSE)
    if(!directed)
      adj <- adj + t(adj)
  }

  ix <- mat2vec.ix(adj, directed, selfloops)
  m <- sum(adj[ix])
  xiregular <- matrix(m^2/sum(adj[ix]!=0), nrow(adj), ncol(adj))
  # if(!directed) xiregular <- xiregular + t(xiregular) - diag(diag(xiregular))
  xiregular <- ceiling(xiregular)

  xiconfiguration <- ComputeXi(adj, directed, selfloops)

  # if(nrow(adj)<50){
    loglikeregular <- extraDistr::dmvhyper(x = adj[ix], n = xiregular[ix], k = m, log = TRUE)
    loglikeconf <- extraDistr::dmvhyper(x = adj[ix], n = xiconfiguration[ix], k = m, log = TRUE)
  # } else{
  #   print(sum(xiregular[ix]))
  #   loglikeregular <- stats::dmultinom(x = adj[ix], prob = xiregular[ix]/sum(xiregular[ix]), log = TRUE)
  #   print(sum(xiconfiguration[ix]))
  #   loglikeconf <- stats::dmultinom(x = adj[ix], prob = xiconfiguration[ix]/sum(xiconfiguration[ix]), log = TRUE)
  # }

  llratio <- -2* (loglikeregular-loglikeconf)

  if(is.null(nempirical)) nempirical <- 100
  ncores <- 1
  if(is.numeric(parallel)) ncores <- parallel
  if(isTRUE(parallel)) ncores <- parallel::detectCores() - 1

  gees <- NULL
  if(nrow(adj)<200){
    gees <- extraDistr::rmvhyper(nn = nempirical, n = xiregular[ix], k = m)
  } else{
    gees <- stats::rmultinom(n = nempirical, prob = xiregular[ix]/sum(xiregular[ix]), size = m)
  }

  nullllratio <- unlist(parallel::mclapply(X = 1:nempirical, FUN = function(id){
    adjr <- NULL
    if(nrow(adj)<200){
      adjr <- gees[id,]
    } else{
      adjr <- gees[,id]
    }
    n <- c(1, nrow(adj),ncol(adj))
    tmp <- vec2mat(adjr,directed,selfloops,n = n)
    if(!directed & n[2]==n[3]) tmp <- tmp + t(tmp) - diag(diag(tmp))
    xiconfigurationr <- ComputeXi(tmp, directed, selfloops)
    loglikeregularr <- extraDistr::dmvhyper(x = adjr, n = xiregular[ix], k = m, log = TRUE)
    loglikeconfr <- extraDistr::dmvhyper(x = adjr, n = xiconfigurationr[ix], k = m, log = TRUE)
    -2 * (loglikeregularr-loglikeconfr)
  }, mc.cores = ncores))

  mm <- -2*m*log(1/sum(ix))
  mu <- mean(nullllratio, na.rm = TRUE)
  va <- var(nullllratio, na.rm = TRUE)

  a <- (mu/(mm*va))*(mu*(mm-mu)-va)
  b <- (mm-mu)*a/mu

  p.value <- pbeta(q = llratio/mm, shape1 = a, shape2 = b, lower.tail = F)
  names(llratio) <- 'lr'
  parms <- nrow(adj)*(1+directed)-1
  names(parms) <- 'df'
  conf.int <- qbeta(p = c(0.025, 0.975), shape1 = a, shape2 = b)*mm
  attributes(conf.int) <- list(conf.level=0.95)
  alternative <- 'one.sided'
  method <- 'LR test -- gnp vs CM'

  return(
    lrtohtest(
      statistic=llratio, parameter=parms, p.value=p.value, conf.int=conf.int,
      alternative=alternative, method=method, data.name=NULL
    )
    )

}


#' Perform likelihood ratio test between two ghype models.
#'
#' lr.test allows to test between two nested ghype models whether there is
#' enough evidence for the alternative (more complex) model compared to the null model.
#'
#' @param nullmodel ghype object. The null model
#' @param altmodel ghype object. The alternative model
#' @param df optional scalar. the number of degrees of freedom.
#' @param williams (deprecated keep FALSE)
#' @param Beta boolean, whether to use empirical Beta distribution approximation. Default TRUE
#' @param seed scalar, seed for the empirical distribution.
#' @param nempirical optional scalar, number of replicates for empirical beta distribution.
#' @param parallel optional, number of cores to use or boolean for parallel computation.
#' If passed TRUE uses all cores-1, else uses the number of cores passed. If none passed
#' performed not in parallel.
#' @param returnBeta boolean, return estimated parameters of Beta distribution? Default FALSE.
#'
#' @return
#'  p-value of test. If returnBeta=TRUE returns the p-value together with the parameters
#'  of the beta distribution.
#' @export
#'
lr.test <- function(nullmodel, altmodel, df=NULL, williams = FALSE, Beta = TRUE, seed = NULL, nempirical = NULL, parallel = FALSE, returnBeta = FALSE){
  llratio <- loglratio(nullmodel,altmodel)

  if(!is.null(Beta)){

    if(is.numeric(Beta)){
      a <- Beta[1]
      b <- Beta[2]
      mm <- Beta[3]

      p.value <- pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F)
      names(llratio) <- 'lr'
      parms <- altmodel$df-nullmodel$df
      names(parms) <- 'df'
      conf.int <- qbeta(p = c(0.025, 0.975), shape1 = a, shape2 = b)*mm
      attributes(conf.int) <- list(conf.level=0.95)
      alternative <- 'one.sided'
      method <- 'LR test'

      return(
        lrtohtest(
          statistic=-2*llratio, parameter=parms, p.value=p.value, conf.int=conf.int,
          alternative=alternative, method=method, data.name=NULL
        )
      )
    }

    if(isTRUE(Beta)){
      if(is.null(nempirical)) nempirical <- 100

      directed <- nullmodel$directed
      selfloops <- nullmodel$selfloops

      ix <- as.matrix(mat2vec.ix(nullmodel$xi,TRUE,TRUE))
      if(length(nullmodel$n)==1)
        ix <- as.matrix(mat2vec.ix(nullmodel$xi,directed,selfloops))
      ps <- nullmodel$omega[ix]*nullmodel$xi[ix]
      mm <- -2*nullmodel$m*log(min(ps[ps!=0])/sum(ps))

      if(is.numeric(seed)) set.seed(seed)

      gees <- rghype(nempirical, nullmodel)

      ncores <- 1
      if(is.numeric(parallel)) ncores <- parallel
      if(isTRUE(parallel)) ncores <- parallel::detectCores() - 1

      nullllratio <- unlist(parallel::mclapply(X = gees, FUN = function(g, directed, null, alt, bip){
        if(!directed &!bip){ # && length(nullmodel$n)==1){
          adj <- g + t(g)
        } else{
          adj <- g
        }
        empnull <- eval(updateModel(nullmodel,adj))
        empalt <- eval(updateModel(altmodel,adj))
        return(-2*loglratio(empnull,empalt))
      }, mc.cores = ncores, directed = directed, null = nullmodel, alt = altmodel, bip = length(nullmodel$n)>1))

      mu <- mean(nullllratio)
      va <- var(nullllratio)

      a <- (mu/(mm*va))*(mu*(mm-mu)-va)
      b <- (mm-mu)*a/mu

      if(returnBeta)
        return(c(pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F), a,b,mm))

      p.value <- pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F)
      names(llratio) <- 'lr'
      parms <- altmodel$df-nullmodel$df
      names(parms) <- 'df'
      conf.int <- qbeta(p = c(0.025, 0.975), shape1 = a, shape2 = b)*mm
      attributes(conf.int) <- list(conf.level=0.95)
      alternative <- 'one.sided'
      method <- 'LR test'

      return(
        lrtohtest(
          statistic=-2*llratio, parameter=parms, p.value=p.value, conf.int=conf.int,
          alternative=alternative, method=method, data.name=NULL
        )
      )
    }
  }

  if(is.null(df)){
    df <- altmodel$df-nullmodel$df
  }

  q1 <- 1
  if(williams){
    ix <- as.matrix(mat2vec.ix(adj,directed,F))
    ps <- nullmodel$omega[ix]*nullmodel$xi[ix]
    ps <- ps[ps!=0]
    k <- df
    q1 <- 1 + ( 6 * nullmodel$m * (k-1) )^(-1) * ( sum( (ps/sum(ps))^(-1) ) - 1 )
  }

  return(pchisq(q = -2*llratio/q1, df = df, lower.tail = FALSE))
}

#' Test SCM vs full ghype.
#'
#' isNetwork tests a graph for the SCM vs the full ghype model.
#'
#' @param graph adjacency matrix or igraph graph
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param nempirical optional, number of graphs to sample from null distribution for empirical distribution.
#' @param parallel optional, number of cores to use or boolean for parallel computation.
#' If passed TRUE uses all cores-1, else uses the number of cores passed. If none passed
#' performed not in parallel.
#' @param returnBeta boolean, return estimated parameters of Beta distribution? Default FALSE.
#'
#'
#' @return
#' p-value of test.
#'
#' @export
isNetwork <- function(graph, directed, selfloops, Beta=NULL, nempirical=NULL, parallel = FALSE, returnBeta = FALSE){
  conftest <- conf.test(graph, directed = directed, selfloops = selfloops, nempirical = nempirical, parallel = parallel)
  if(conftest$p.value >= 1e-3){
    adj <- graph
    if(requireNamespace("igraph", quietly = TRUE) && igraph::is.igraph(graph)){
      adj <- igraph::get.adjacency(graph, type='upper', sparse = FALSE)
      if(!directed)
        adj <- adj + t(adj)
    }
    ix <- mat2vec.ix(adj, directed, selfloops)
    m <- sum(adj[ix])
    xiregular <- matrix(m^2/sum(adj[ix]!=0), nrow(adj), ncol(adj))
    # if(!directed) xiregular <- xiregular + t(xiregular) - diag(diag(xiregular))
    xiregular <- ceiling(xiregular)
    fullmod <- ghype(graph, directed, selfloops)
    nullmod <- ghype(object = graph, directed = directed, selfloops = selfloops, xi = xiregular, unbiased = TRUE)
  } else{
    fullmod <- ghype(graph, directed, selfloops)
    nullmod <- ghype(graph, directed, selfloops, unbiased = TRUE)
  }
  # n <- full$n[1]
  # df <- n*(n-!selfloops)/(1+!directed)
  # if(igraph::is.igraph(graph)){
  #   if(igraph::is.bipartite(graph))
  #     df <- sum(V(graph)$type)*sum(!V(graph)$type)
  # }
  #
  # if(is.matrix(graph)){
  #   if(nrow(graph)!=ncol(graph)){
  #     df <- nrow(graph)*ncol(graph)
  #   }
  # }
  return(lr.test(nullmodel = nullmod,altmodel = fullmod, Beta=Beta, nempirical = nempirical, parallel = parallel, returnBeta = returnBeta))
}


#' Estimate statistical deviations from ghype model
#'
#' linkSignificance allows to estimate the statistical deviations of an observed
#' graph from a ghype model.
#'
#' @param graph an adjacency matrix or a igraph object.
#' @param model a ghype model
#' @param under boolean, estimate under-represented deviations? Default FALSE.
#'
#' @return
#'
#' matrix of probabilities with same size as adjacency matrix.
#'
#' @export
#'
linkSignificance <- function(graph, model, under=FALSE){
  adj <- graph
  if(requireNamespace("igraph", quietly = TRUE) && igraph::is.igraph(graph)){
    adj <- igraph::get.adjacency(graph, type='upper', sparse = FALSE)
    if(!directed)
      adj <- adj + t(adj)
  }
  graph <- adj

  directed <- model$directed
  selfloops <- model$selfloops

  # get relevant indices
  idx <- mat2vec.ix(graph, directed, selfloops)

  # compute parameters for marginal distributions
  xibar <- sum(model$xi[idx])-model$xi[idx]
  omegabar <- (sum(model$xi[idx]*model$omega[idx])-model$xi[idx]*model$omega[idx])/xibar

  # compute vector of probabilities using Wallenius univariate distribution or binomial
  id <- graph[idx]!=0
  probvec <- rep(1, sum(idx))

  if((mean(xibar)/sum(graph[idx]))<1e3){
    probvec[id] <- Vectorize(FUN = BiasedUrn::pWNCHypergeo, vectorize.args = c('x', 'm1', 'm2','n','odds'))(
      x = graph[idx][id],m1 = model$xi[idx][id],m2 = xibar[id],
      n = sum(graph[idx]), odds = model$omega[idx][id]/omegabar[id],
      lower.tail = under
      )
  } else{
    probvec[id] <- Vectorize(FUN = stats::pbinom, vectorize.args = c('q', 'size', 'prob'))(
      q = graph[idx][id], size = sum(graph[idx]),
      prob = model$xi[idx][id]* model$omega[idx][id]/(
            model$xi[idx][id] * model$omega[idx][id]+xibar[id]*omegabar[id]
            ),
      lower.tail = under
      )
  }

  # return matrix of significance for each entry of original adjacency
  return(vec2mat(probvec,directed,selfloops,nrow(graph)))
}
