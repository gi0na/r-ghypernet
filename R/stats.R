# Computes Fisher Information matrix for estimators in nrm models.
Jn <- function(beta=NULL, w=NULL, xi=NULL, adj=NULL, 
    directed=NULL, selfloops=NULL, model=NULL) {
    # Returns Fisher Information
    # matrix w has to be a list to
    # be passed to an apply function
  
    if(!is.null(model)){
      beta <- coef(model)
      xi <- model$xi
      directed <- model$directed
      selfloops <- model$selfloops
      w <- model$predictors
      adj <- model$adj
    }
  
    if (!is.list(w)) 
        w <- list(w)
    # take upper triangular part of
    # adjacency matrix
    ix <- mat2vec.ix(adj, directed, 
        selfloops)
    # compute the product of all
    # layers to the power of the
    # parameter beta
    product <- apply(sapply(X = 1:length(beta), 
        FUN = function(i) {
            w[[i]][ix]^beta[i]
        }), MARGIN = 1, FUN = prod)
    # return the value of the rhs in
    # eq. 7
    a <- sum(xi[ix] * product)
    
    b <- sapply(X=1:length(beta),
                FUN=function(i){
                  sapply(X=1:length(beta), FUN = function(j,i){
                    sum(xi[ix]*product*log(w[[i]][ix])*log(w[[j]][ix]))
                  }, i=i)
                })
    
    c <- sapply(X=1:length(beta),
               FUN=function(i){
                 sum(log(w[[i]][ix])*xi[ix]*product)
               })
    
    d <- tcrossprod(c)
    
    return(sum(adj[ix], na.rm = TRUE) * 
        (a*b-d)/(a^2))
}


#' Computes Mc Fadden pseudo R-squared.
#' 
#' Pass either the models or the model parameters as arguments
#' 
#' @param adj optimal adjacency matrix
#' @param xi  optional xi matrix
#' @param omega0  optional propensity matrix of null model
#' @param omega1  optional propensity matrix of alternative model
#' @param directed  boolean, is the model directed?
#' @param selfloops  boolean, are there selfloops?
#' @param mod0  nrm null model
#' @param mod1  nrm alternative model
#' @param nparam  integer, number of parameters
#' @return Mc Fadden pseudo R-squared.
#' 
#' @export
mcfaddenR2 <- function(adj = NULL, 
                       xi = NULL, omega0 = NULL, omega1 = NULL, 
                       directed, selfloops, 
                       mod0 = NULL, mod1 = NULL, nparam) {
  # Returns McFadden adjusted
  # pseudo-R-squared (McFadden,
  # 1974).
  if (is.null(mod0) | is.null(mod1)) {
    R2 <- 1 - ((logl(object = adj, 
                     xi, omega = omega1, 
                     directed = directed, selfloops = selfloops) - 
                  nparam)/logl(object = adj, 
                               xi = xi, omega = omega0, 
                               directed = directed, selfloops = selfloops))
  } else {
    R2 <- 1 - ((mod1$loglikelihood - 
                  nparam)/mod0$loglikelihood)
  }
  return(R2)
}


#' Computes Cox and Snell pseudo R-squared for nrm models.
#' 
#' @param mod0  nrm null model
#' @param mod1  nrm alternative model
#' @param m  number of edges
#' @return Cox and Snell pseudo R-squared
#' @author  GC
#' 
#' @export
coxsnellR2 <- function(mod0, mod1, m) {
  # Returns Cox-Snell
  # pseudo-R-squared.
  -expm1(2 * (mod0$loglikelihood - 
                mod1$loglikelihood)/m)
}

# vcov
#' @exportS3Method
vcov.nrm <- function(object, ...){
  solve(Jn(model = object))
}
