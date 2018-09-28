#' Title
#'
#' @param graph
#' @param directed
#' @param selfloops
#' @param nempirical
#'
#' @return
#' @export
#'
#' @examples
configurationtest <- function(graph, directed, selfloops, nempirical=NULL){

  adj <- graph
  if(requireNamespace("igraph", quietly = TRUE) && igraph::is.igraph(graph)){
    adj <- igraph::get.adjacency(graph, type='upper', sparse = FALSE)
    if(!directed)
      adj <- adj + t(adj)
  }

  ix <- mat2vec.ix(adj, directed, selfloops)
  m <- sum(adj[ix])
  xiregular <- matrix(m^2/sum(ix), nrow(adj), ncol(adj))
  if(!directed) xiregular <- xiregular + t(xiregular) - diag(diag(xiregular))
  xiregular <- round(xiregular)

  xiconfiguration <- ComputeXi(adj, directed, selfloops)

  loglikeregular <- extraDistr::dmvhyper(x = adj[ix], n = xiregular[ix], k = m, log = TRUE)
  loglikeconf <- extraDistr::dmvhyper(x = adj[ix], n = xiconfiguration[ix], k = m, log = TRUE)
  llratio <- -2* (loglikeregular-loglikeconf)

  if(is.null(nempirical)) nempirical <- 100

  gees <- NULL
  if(nrow(adj)<200){
    gees <- extraDistr::rmvhyper(nn = nempirical, n = xiregular[ix], k = m)
  } else{
    gees <- stats::rmultinom(n = nempirical, prob = xiregular[ix]/sum(xiregular[ix]), size = m)
  }

  nullllratio <- sapply(X = 1:nempirical, FUN = function(id){
    adjr <- NULL
    if(nrow(adj)<200){
      adjr <- gees[id,]
    } else{
      adjr <- gees[,id]
    }
    tmp <- vec2mat(adjr,directed,selfloops,n = nrow(adj))
    if(!directed) tmp <- tmp + t(tmp) - diag(diag(tmp))
    xiconfigurationr <- ComputeXi(tmp, directed, selfloops)
    loglikeregularr <- extraDistr::dmvhyper(x = adjr, n = xiregular[ix], k = m, log = TRUE)
    loglikeconfr <- extraDistr::dmvhyper(x = adjr, n = xiconfigurationr[ix], k = m, log = TRUE)
    -2 * (loglikeregularr-loglikeconfr)
  })

  mm <- -2*m*log(1/sum(ix))
  mu <- mean(nullllratio)
  va <- var(nullllratio)

  a <- (mu/(mm*va))*(mu*(mm-mu)-va)
  b <- (mm-mu)*a/mu

  return(pbeta(q = llratio/mm, shape1 = a, shape2 = b, lower.tail = F))

}


#' Title
#'
#' @param nullmodel
#' @param altmodel
#' @param df
#'
#' @return
#' @export
#'
#' @examples
llratiotest <- function(nullmodel, altmodel, df=NULL, williams = FALSE, Beta = TRUE, seed = NULL, nempirical = NULL, parallel = FALSE, returnBeta = FALSE){
  llratio <- loglratio(nullmodel,altmodel)

  if(!is.null(Beta)){

    if(is.numeric(Beta)){
      a <- Beta[1]
      b <- Beta[2]
      mm <- Beta[3]
      return(pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F))
    }

    if(isTRUE(Beta)){
      if(is.null(nempirical)) nempirical <- 100

      directed <- nullmodel$directed
      selfloops <- nullmodel$selfloops

      ix <- as.matrix(mat2vec.ix(nullmodel$xi,directed,selfloops))
      ps <- nullmodel$omega[ix]*nullmodel$xi[ix]
      mm <- -2*nullmodel$m*log(min(ps[ps!=0])/sum(ps))

      if(is.numeric(seed)) set.seed(seed)

      gees <- RandomGraph(nempirical, nullmodel)

      ncores <- 1
      if(is.numeric(parallel)) ncores <- parallel
      if(isTRUE(parallel)) ncores <- parallel::detectCores() - 1

      nullllratio <- unlist(parallel::mclapply(X = gees, FUN = function(g, directed, null, alt){
        if(!directed){
          adj <- g + t(g)
        } else{
          adj <- g
        }
        empnull <- eval(updateModel(null,adj))
        empalt <- eval(updateModel(alt,adj))
        return(-2*loglratio(empnull,empalt))
      }, mc.cores = ncores, directed = directed, null = nullmodel, alt = altmodel))

      mu <- mean(nullllratio)
      va <- var(nullllratio)

      a <- (mu/(mm*va))*(mu*(mm-mu)-va)
      b <- (mm-mu)*a/mu

      if(returnBeta)
        return(c(pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F), a,b,mm))
      return(pbeta(q = -2*llratio/mm, shape1 = a, shape2 = b, lower.tail = F))
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

#' Title
#'
#' @param graph
#' @param directed
#' @param selfloops
#'
#' @return
#' @export
#'
#' @examples
isNetwork <- function(graph, directed, selfloops, Beta=NULL, nempirical=NULL, parallel = FALSE, returnBeta = FALSE){
  fullmod <- ghype(graph, directed, selfloops)
  nullmod <- ghype(graph, directed, selfloops, unbiased = TRUE)
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
  return(llratiotest(nullmodel = nullmod,altmodel = fullmod, Beta=Beta, nempirical = nempirical, parallel = parallel, returnBeta = returnBeta))
}

#' Title
#'
#' @param graph
#' @param directed
#' @param selfloops
#'
#' @return
#' @export
#'
#' @examples
linkSignificance <- function(graph, model){
  ## TODO: assume graph is adjacency for now. Extend with method for graph and for matrices

  directed <- model$directed
  selfloops <- model$selfloops

  # get relevant indices
  idx <- mat2vec.ix(graph, directed, selfloops)

  # compute parameters for marginal distributions
  xibar <- sum(model$xi[idx])-model$xi[idx]
  omegabar <- (sum(model$xi[idx]*model$omega[idx])-model$xi[idx]*model$omega[idx])/xibar

  # compute vector of probabilities using Wallenius univariate distribution
  probvec <- Vectorize(FUN = function(id){
    BiasedUrn::pWNCHypergeo(x = graph[idx][id],m1 = model$xi[idx][id],m2 = xibar[id], n = sum(graph[idx]), odds = model$omega[idx][id]/omegabar[id])
    }, vectorize.args = 'id')(1:sum(idx))

  # return matrix of significance for each entry of original adjacency
  return(vec2mat(probvec,directed,selfloops,nrow(graph)))
}
