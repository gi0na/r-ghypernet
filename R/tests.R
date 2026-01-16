#' Test regular (gnp) vs configuration model
#'
#' Likelihood ratio test for gnp vs configuration model.
#'
#' @importFrom methods new 
#' @param graph adjacency matrix or igraph graph
#' @param directed a boolean argument specifying whether object is directed or not.
#' @param selfloops a boolean argument specifying whether the model should incorporate selfloops.
#' @param nempirical optional, number of graphs to sample from null distribution for empirical distribution.
#' @param parallel optional, number of cores to use or boolean for parallel computation.
#' If passed TRUE uses all cores-1, else uses the number of cores passed. If none passed
#' performed not in parallel.
#' @param seed optional integer
#'
#' @return
#' p-value of test.
#'
#' @export
#'
#' @examples
#' data("adj_karate")
#' conf.test(graph = adj_karate, directed = FALSE, selfloops = FALSE, seed=123)
#' 
conf.test <- function(graph, directed, selfloops, nempirical=NULL, parallel=NULL, seed = NULL){
  if(is.numeric(seed)){
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
    set.seed(seed)
  }

  adj <- graph
  if(requireNamespace("igraph", quietly = TRUE) && igraph::is_igraph(graph)){
    adj <- igraph::get.adjacency(graph, type='upper', sparse = FALSE)
    if(!directed)
      adj <- adj + t(adj)
  }

  ix <- mat2vec.ix(adj, directed, selfloops)
  m <- sum(adj[ix])
  xiregular <- matrix(m^2/sum(adj[ix]!=0), nrow(adj), ncol(adj))
  # if(!directed) xiregular <- xiregular + t(xiregular) - diag(diag(xiregular))
  xiregular <- ceiling(xiregular)

  xiconfiguration <- compute_xi(adj, directed, selfloops)

  # if(nrow(adj)<50){
    loglikeregular <- dmvhyper_base(x = adj[ix], n = xiregular[ix], k = m, log = TRUE)
    loglikeconf <- dmvhyper_base(x = adj[ix], n = xiconfiguration[ix], k = m, log = TRUE)
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
    gees <- rmvhyper_base(nn = nempirical, n = xiregular[ix], k = m)
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
    xiconfigurationr <- compute_xi(tmp, directed, selfloops)
    loglikeregularr <- dmvhyper_base(x = adjr, n = xiregular[ix], k = m, log = TRUE)
    loglikeconfr <- dmvhyper_base(x = adjr, n = xiconfigurationr[ix], k = m, log = TRUE)
    -2 * (loglikeregularr-loglikeconfr)
  }, mc.cores = ncores))

  mm <- -2*m*log(1/sum(ix))
  mu <- mean(nullllratio, na.rm = TRUE)
  va <- stats::var(nullllratio, na.rm = TRUE)

  a <- (mu/(mm*va))*(mu*(mm-mu)-va)
  b <- (mm-mu)*a/mu

  p.value <- stats::pbeta(q = llratio/mm, shape1 = a, shape2 = b, lower.tail = F)
  names(llratio) <- 'lr'
  parms <- nrow(adj)*(1+directed)-1
  names(parms) <- 'df'
  conf.int <- stats::qbeta(p = c(0.025, 0.975), shape1 = a, shape2 = b)*mm
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
# @param williams (deprecated keep FALSE)
#' @param Beta boolean, whether to use empirical Beta distribution approximation. Default TRUE
#' @param seed scalar, seed for the empirical distribution.
#' @param nempirical optional scalar, number of replicates for empirical beta distribution.
#' @param parallel optional, number of cores to use or boolean for parallel computation.
#' If passed TRUE uses all cores-1, else uses the number of cores passed. If none passed
#' performed not in parallel.
#' @param returnBeta boolean, return estimated parameters of Beta distribution? Default FALSE.
#' @param method string, for internal use
#'
#' @return
#'  p-value of test. If returnBeta=TRUE returns the p-value together with the parameters
#'  of the beta distribution.
#' @export
#' 
#' @examples
#' data("adj_karate")
#' regularmodel <- regularm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' confmodel <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' lr.test(nullmodel = regularmodel, altmodel = confmodel, seed = 123)
#'
lr.test <- function(nullmodel, altmodel, df=NULL, Beta = TRUE, seed = NULL, nempirical = NULL, parallel = FALSE, returnBeta = FALSE, method = NULL){
  llratio <- loglratio(nullmodel,altmodel)
  if(is.numeric(seed)){
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
    set.seed(seed)
  }

  if(is.null(method))
    method <- 'LR test'

  if(!is.null(Beta)){

    if(is.numeric(Beta)){
      a <- Beta[1]
      b <- Beta[2]
      mm <- Beta[3]

      p.value <- stats::pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F)
      names(llratio) <- 'lr'
      parms <- altmodel$df-nullmodel$df
      names(parms) <- 'df'
      conf.int <- stats::qbeta(p = c(0.025, 0.975), shape1 = a, shape2 = b)*mm
      attributes(conf.int) <- list(conf.level=0.95)
      alternative <- 'one.sided'

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

      gees <- rghype(nempirical, nullmodel)

      ncores <- 1
      if(is.numeric(parallel)) ncores <- parallel
      if(isTRUE(parallel)) ncores <- parallel::detectCores() - 1

      nullllratio <- unlist(parallel::mclapply(X = gees, FUN = function(adj, directed, null, alt, bip){
        empnull <- eval(updateModel(nullmodel,adj))
        empalt <- eval(updateModel(altmodel,adj))
        return(-2*loglratio(empnull,empalt))
      }, mc.cores = ncores, directed = directed, null = nullmodel, alt = altmodel, bip = length(nullmodel$n)>1))


      mu <- mean(nullllratio)
      va <- stats::var(nullllratio)

      a <- (mu/(mm*va))*(mu*(mm-mu)-va)
      b <- (mm-mu)*a/mu

      if(isTRUE(returnBeta))
        return(c(stats::pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F), a,b,mm))

      p.value <- stats::pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F)
      names(llratio) <- 'lr'
      parms <- altmodel$df-nullmodel$df
      names(parms) <- 'df'
      conf.int <- stats::qbeta(p = c(0.025, 0.975), shape1 = a, shape2 = b)*mm
      attributes(conf.int) <- list(conf.level=0.95)
      alternative <- 'one.sided'

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
  if(FALSE){
    ix <- as.matrix(mat2vec.ix(nullmodel$xi,nullmodel$directed,nullmodel$selfloops))
    ps <- nullmodel$omega[ix]*nullmodel$xi[ix]
    ps <- ps[ps!=0]
    k <- df
    q1 <- 1 + ( 6 * nullmodel$m * (k-1) )^(-1) * ( sum( (ps/sum(ps))^(-1) ) - 1 )
  }

  return(stats::pchisq(q = -2*llratio/q1, df = df, lower.tail = FALSE))
}

#' Test null model vs full ghype.
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
#' @param Beta boolean, use Beta test? default TRUE
#' @param seed optional integer, seed for empirical lr.test
#'
#' @return
#' p-value of test.
#'
#' @export
#' 
#' @examples
#' data("adj_karate")
#' isNetwork(graph = adj_karate, directed = FALSE, selfloops = FALSE, seed=123)
#' 
isNetwork <- function(graph, directed, selfloops, Beta=TRUE, nempirical=NULL, parallel = FALSE, returnBeta = FALSE, seed = NULL){
  conftest <- conf.test(graph, directed = directed, selfloops = selfloops, nempirical = nempirical, parallel = parallel, seed = NULL)
  if(conftest$p.value >= 1e-3){
    method <- 'LR test -- optimal = gnp vs full model'
    adj <- graph
    if(requireNamespace("igraph", quietly = TRUE) && igraph::is_igraph(graph)){
      adj <- igraph::get.adjacency(graph, type='upper', sparse = FALSE)
      if(!directed)
        adj <- adj + t(adj)
    }
    ix <- mat2vec.ix(adj, directed, selfloops)
    m <- sum(adj[ix])
    xiregular <- matrix(m^2/sum(adj[ix]!=0), nrow(adj), ncol(adj))
    # if(!directed) xiregular <- xiregular + t(xiregular) - diag(diag(xiregular))
    xiregular <- ceiling(xiregular)
    fullmod <- ghype(graph, directed, selfloops, xi = xiregular)
    nullmod <- ghype(graph = graph, directed = directed, selfloops = selfloops, xi = xiregular, unbiased = TRUE)
    nullmod$df <- 1
  } else{
    method <- 'LR test -- optimal = CM vs full model'
    fullmod <- ghype(graph, directed, selfloops)
    nullmod <- ghype(graph, directed, selfloops, unbiased = TRUE)
    nullmod$df <- nullmod$n * (1+directed)
  }
  return(lr.test(nullmodel = nullmod,altmodel = fullmod, Beta=Beta, nempirical = nempirical, parallel = parallel, returnBeta = returnBeta, method = method, seed = seed))
}

#' Perform a goodness-of-fit test
#'
#' @param model ghype model to test
#' @param Beta boolean, whether to use empirical Beta distribution approximation. Default TRUE
#' @param nempirical optional scalar, number of replicates for empirical beta distribution.
#' @param parallel optional, number of cores to use or boolean for parallel computation.
#' If passed TRUE uses all cores-1, else uses the number of cores passed. If none passed
#' performed not in parallel.
#' @param returnBeta boolean, return estimated parameters of Beta distribution? Default FALSE.
#' @param seed scalar, seed for the empirical distribution.
#'
#' @return
#'  p-value of test. If returnBeta=TRUE returns the p-value together with the parameters
#'  of the beta distribution.
#' @export
#' 
#' @examples
#' data("adj_karate")
#' confmodel <- scm(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' gof.test(model = confmodel, seed = 123)
#'
gof.test <- function(model, Beta=TRUE, nempirical = NULL, parallel = NULL, returnBeta = FALSE, seed = NULL){
  fullmodel <- ghype(graph = model$adj, directed = model$directed, selfloops = model$selfloops, unbiased = FALSE)
  return(lr.test(nullmodel = model,altmodel = fullmodel, Beta=Beta, nempirical = nempirical, parallel = parallel, returnBeta = returnBeta, seed = seed, method = 'LR test -- GOF'))
}


#' Estimate statistical deviations from ghype model
#'
#' linkSignificance allows to estimate the statistical deviations of an observed
#' graph from a ghype model.
#'
#' @param graph an adjacency matrix or a igraph object.
#' @param model a ghype model
#' @param under boolean, estimate under-represented deviations? Default FALSE:
#'   i.e. returns over representation
#' @param log.p boolean, return log values of probabilities
#' @param binomial.approximation boolean, force binomial? default FALSE
#' @param give_pvals boolean, return p-values for both under and over
#'   significance? when FALSE, it returns probabilty of observing stricly more
#'   (or less) edges than in graph. When TRUE returns probability of observing
#'   exactly as many edges or more (less) than in graph, like a standard pvalue.
#'
#' @return
#'
#' matrix of probabilities with same size as adjacency matrix.
#'
#' @export
#'
#' @examples
#' data("adj_karate")
#' fullmodel <- ghype(graph = adj_karate, directed = FALSE, selfloops = FALSE)
#' link_significance(graph = adj_karate, model = fullmodel, under=FALSE)
#' 
link_significance <- function(graph, model, under=FALSE, log.p=FALSE, binomial.approximation = FALSE, give_pvals = TRUE){
  adj <- graph
  if(requireNamespace("igraph", quietly = TRUE) && igraph::is_igraph(graph)){
    adj <- igraph::get.adjacency(graph, type='both', sparse = FALSE)
    # if(!directed)
    #   adj <- adj + t(adj)
  }

  directed <- model$directed
  selfloops <- model$selfloops

  # get relevant indices
  idx <- mat2vec.ix(adj, directed, selfloops)

  # compute parameters for marginal distributions
  xibar <- sum(model$xi[idx])-model$xi[idx]
  omegabar <- (sum(model$xi[idx]*model$omega[idx])-model$xi[idx]*model$omega[idx])/xibar

  # compute vector of probabilities using hypergeometric, Wallenius univariate distribution or binomial
  if(!under & give_pvals){
    id <- adj[idx]!=0
  } else{
    id <- is.numeric(adj[idx])
  }
  probvec <- rep(ifelse(log.p, 0, 1), sum(idx))

  if( all(model$omega[idx] == model$omega[1]) & (!binomial.approximation) ){
    probvec[id] <- stats::phyper(
      q = adj[idx][id] - (give_pvals&!under) - (!give_pvals&under), m = model$xi[idx][id], n = xibar[id],
      k = sum(adj[idx]),
      lower.tail = under, log.p = log.p
    )
  } else{

    if( requireNamespace("BiasedUrn", quietly = TRUE) && (((mean(xibar)/sum(adj[idx]))<1e3) & (!binomial.approximation)) ){
      probvec[id] <- mapply(FUN = BiasedUrn::pWNCHypergeo,
        x = adj[idx][id]  - (give_pvals&!under) - (!give_pvals&under), m1 = model$xi[idx][id],m2 = xibar[id],
        odds = model$omega[idx][id]/omegabar[id], MoreArgs = list(
        n = sum(adj[idx]), lower.tail = under
        ))
      if(log.p) probvec[id] <- log(probvec[id])
    } else{
      probvec[id] <- stats::pbinom(
        q = adj[idx][id]  - (give_pvals&!under) - (!give_pvals&under), size = sum(adj[idx]),
        prob = model$xi[idx][id]* model$omega[idx][id]/(
              model$xi[idx][id] * model$omega[idx][id]+xibar[id]*omegabar[id]
              ),
        lower.tail = under, log.p = log.p
        )
    }
  }

  # return matrix of significance for each entry of original adjacency
  return(vec2mat(probvec, directed, selfloops, model$n))
}


