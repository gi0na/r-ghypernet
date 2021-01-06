updateModel <- function(model, adj){

  callname <- as.character(model$call[1])

  xi <- NULL

  newcall <- NULL

  if(length(grep('function', callname))>0){
    if(length(grep('block', callname))>0 | length(grep('labels', callname))>0){
      callname <- 'bccm'
    } else{
      callname <- 'ghype'
    }
  } else{
    fixXi <- length(grep('xi', deparse(model$call, width.cutoff = 500)))>0
    if(fixXi)
      xi <- model$xi
  }
  if(length(grep('ghype', callname))>0){
    callname <- 'ghype'
    newcall <- call(name = callname, graph=adj, directed=model$directed, selfloops=model$selfloops, xi=xi, unbiased=all(model$omega==1), regular=model$regular, multinomial = TRUE)
  } else{
    if(length(grep('bccm', callname))>0){
      callname <- 'bccm'
      newcall <- call(name = callname, adj=adj, labels=model$labels, directed=model$directed, selfloops=model$selfloops, xi=xi, regular=model$regular, directedBlocks=model$directedBlocks, homophily=model$homophily,inBlockOnly=model$inBlockOnly, multinomial = TRUE)
    }
    if(length(grep('nrm', callname))>0){
      callname <- 'nrm'
      newcall <- call(name = callname, w=model$predictors, adj=adj, directed=model$directed,
                      selfloops=model$selfloops, xi=xi, init = c(model$coef[-length(model$coef)],0.01), ci=FALSE, regular=model$regular, multinomial = TRUE)
    }
  }

  return(newcall)
}

lrtohtest <- function(statistic, parameter, p.value, conf.int, alternative, method, data.name){
  val <- list(statistic=statistic, parameter=parameter,
              p.value=p.value, conf.int=conf.int, alternative=alternative,
              method=method, data.name=data.name)
  class(val) <- 'htest'
  return(val)
}

check_specs <- function(object, ...){
  UseMethod('check_specs')
}

check_specs.matrix <- function(object, ...){
  if(is.matrix(object)){
    # if(is.null(directed)){
      if(isSymmetric(object)){
        directed <- FALSE
      } else{
        directed <- TRUE
      }
    # } else{
    #    if(!directed & !isSymmetric(object)){
    #     warning('Trying to compute undirected ensemble for asymmetric adjacency matrix.
    #           Adjacency matrix symmetrised as adj <- adj + t(adj)')
    #     object <- object + t(object)
    #   }
    # }

    # if(is.null(selfloops)){
      if(all(diag(object)==0)){
        selfloops <- FALSE
      } else{
        selfloops <- TRUE
      }
    # }
  }
  return(c('directed'=directed, 'selfloops'=selfloops))
}