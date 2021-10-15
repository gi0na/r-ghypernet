#' Fitting gHypEG regression models for multi-edge networks.
#'
#' nrm is used to fit multi-edge network regression models.
#'
#' @param w an object of class \code{'list'} containing the predictors layers
#'   (explanatory variables/covariates) of the multiplex, passed as adjacency
#'   matrices. The entries of the list can be named.
#' @param adj matrix. The adjacency matrix of the response network (dependent
#'   variable).
#' @param xi optional matrix. Passes a non-standard \eqn{\Xi} matrix.
#' @param pval the significance level used to compute confidence intervals of
#'   the parameters. Per default, set to 0.01.
#' @param directed logical. If \code{TRUE} the response variable is considered
#'   the adjacency matrix of directed graph.  If \code{FALSE} only the upper
#'   triangular of \code{adj} is considered. Default set to FALSE.
#' @param selfloops logical. Whether selfloops are allowed. Default set to
#'   FALSE.
#' @param regular logical. Whether the gHypEG regression should be performed
#'   with correction of combinatorial effects (\code{TRUE}) or without
#'   (\code{FALSE}).
#' @param ci logical. Whether to compute confidences for the parameters.
#'   Defaults to \code{TRUE}.
#' @param significance logical. Whether to test the model significance against
#'   the null by means of lr-test.
#' @param null logical. Is this a null model? Used for internal routines.
#' @param init numeric. Vector of initial values used for numerical MLE. If only
#'   a single value is passed, this is repeated to match the number of
#'   predictors in \code{w}.
#' @param \dots additional arguments to be passed to the low level regression
#'   fitting functions.
#' @return nrm returns an object of class 'nrm'.
#'
#'   The function summary is used to obtain and print a summary and analysis of
#'   the results. The generic accessory functions coefficients, etc, extract
#'   various useful features of the value returned by nrm.
#'
#'   An object of class 'nrm' is a list containing at least the following
#'   components:
#'
#'   \item{coef }{a named vector of coefficients.} \item{confint }{a named
#'   matrix with confidence intervals and standard deviation for each
#'   coefficient.} \item{omega }{the estimated propensity matrix.} \item{xi
#'   }{the matrix of possibilities.} \item{loglikelihood }{log-likelihood of the
#'   estimated model.} \item{AIC }{AIC of the estimated model.} \item{R2 }{Mc
#'   Fadden pseudo R-squared} \item{csR2 }{Cox and Snells pseudo R-squared}
#'   \item{significance }{the p-value of the likelihood-ratio test for the
#'   estimated model against the null.}
#' @author Giona Casiraghi
#' @references Casiraghi, Giona. 'Multiplex Network Regression: How do relations
#'   drive interactions?.' arXiv preprint arXiv:1702.02048 (2017).
#' @keywords models regression nonlinear multivariate sna
#' @examples
#'
#' ## For a complete example see the vignette
#' 
#' data('highschool.predictors')
#'
#' highschool.m <- nrm(w=highschool.predictors[1], adj=contacts.adj, directed=FALSE,
#'   selfloops=FALSE)
#'
#' highschool.m
#' 
#' \donttest{
#' data('highschool.predictors')
#'
#' highschool.m <- nrm(w=highschool.predictors, adj=contacts.adj, directed=FALSE,
#'   selfloops=FALSE)
#'
#' highschool.m
#'}
#'
#' @export
nrm <- function(w, adj, xi = NULL, 
    pval = 0.01, directed = TRUE, 
    selfloops = TRUE, regular = FALSE,
    ...) UseMethod("nrm")

#' @describeIn nrm Default method for nrm
#' @export
#'
nrm.default <- function(w, adj, 
                        xi = NULL, pval = 0.01, directed = FALSE, 
                        selfloops = FALSE, regular = FALSE, ci = TRUE, significance = FALSE, 
                        null = FALSE, init = NULL, ...) {
  # Estimate the multivariate
  # network regression model
  # returns a list with
  # coefficients, joint omega
  # matrix, log-likelihood value,
  # AIC, and R2
  if (is.null(xi)) {
    ## Build Ensemble
    xi <- compute_xi(adj, directed, 
                    selfloops, regular)
  }
  if (!null) {
    mod0 <- nrm.default(w = list(matrix(1, 
                                        nrow(adj), ncol(adj))), 
                        adj = adj, xi, directed = directed, 
                        selfloops = selfloops, 
                        ci = FALSE, significance = FALSE, 
                        null = TRUE)
  }
  # number parameters for
  # computing the xi matrix
  k.xi <- nrow(xi)*(1+directed)*(1-regular) + regular
  # MLE of beta parameters
  if (length(w) == 1 && length(unique(as.vector(w[[1]]))) == 
      1) {
    b <- 0
  } else {
    if (is.null(init)) {
      init <- rep(0.9, length(w))
    } else {
      init <- c(init, rep(0.9, 
                          length(w) - length(init)))
    }
    b <- rootSolve::multiroot(f = fnM, 
                              jacfunc = Jn, w = w, 
                              xi = xi, adj = adj, 
                              directed = directed, 
                              selfloops = selfloops, 
                              start = unlist(init), rtol = 1e-10)$root
  }
  # omega matrix
  omega <- apply(sapply(X = 1:length(b), 
                        FUN = function(i) {
                          w[[i]]^b[i]
                        }), MARGIN = 1, FUN = prod)
  # Compute the various stats
  ll <- logl(object = adj, xi = xi, 
             omega = omega, directed=directed, 
             selfloops = selfloops)
  # AIC <-
  # (log(sum(mat2vec.ix(adj,directed,selfloops)))
  # * (length(b) + k.xi) - 2 * ll)
  AIC <- (2 * (length(b) + k.xi) - 
            2 * ll)
  DL <- -ll/log(2) + (length(b) + k.xi)/2*
    log(sum(adj[mat2vec.ix(adj,directed,selfloops)]), base = 2)
  R2 <- mcfaddenR2(adj = adj, 
                   xi = xi, omega0 = matrix(1, 
                                            nrow(adj), ncol(adj)), 
                   omega1 = omega, directed = directed, 
                   selfloops = selfloops, nparam = length(b) + 
                     k.xi)
  names(b) <- names(w)
  
  # TODO: build from ghype
  mod <- list(call = match.call(), 
              coef = b, confint = NULL, 
              omega = omega, xi = xi, 
              loglikelihood = ll, AIC = AIC, DL = DL,
              R2 = R2, csR2 = 0, directed = directed, 
              selfloops = selfloops, pvalue = pval, 
              significance = NULL,
              adj = adj, predictors = w, df=length(b) + k.xi,
              m = sum(adj[mat2vec.ix(adj, directed, selfloops)]), n = nrow(xi),
              regular = regular)
  if (!null) {
    mod$csR2 <- coxsnellR2(mod0 = mod0, 
                           mod1 = mod, 
                           m = sum(adj[mat2vec.ix(adj, 
                                                       directed, selfloops)]))
  }
  if (ci) {
    mod$confint <- nr.ci(nr.m = mod, 
                         w = w, adj = adj, pval = pval)
  } else{
    mod$confint <- matrix(NA,1,3)
  }
  if (significance) {
    mod$significance <- nr.significance(mod1 = mod, 
                                        adj = adj)
  }
  class(mod) <- c("nrm", 'ghype')
  return(mod)
}


# ' Auxilliary function. Returns a value proportional to the first derivative of
# ' the likelihood in nrm models.
# ' 
# '  ~~ A concise (1-5 lines) description of what the function does. ~~
# ' 
# '  ~~ If necessary, more details than the description above ~~
# ' 
# ' @param x  ~~Describe \code{x} here~~
# ' @param w  ~~Describe \code{w} here~~
# ' @param xi  ~~Describe \code{xi} here~~
# ' @param adj  ~~Describe \code{adj} here~~
# ' @param directed  ~~Describe \code{directed} here~~
# ' @param selfloops  ~~Describe \code{selfloops} here~~
# ' @return val
# ' @note  ~~further notes~~
# ' @author  ~~who you are~~
# ' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
# ' @references  ~put references to the literature/web site here ~
# ' @keywords ~kwd1 ~kwd2
# ' @examples
fnM <- function(x, w, xi, adj, directed, 
                selfloops) {
  # Returns the term of the
  # negative log-likelihood to be
  # minimised (cf.  eq.7).  It
  # assumes the network is
  # undirected and without
  # selfloops.  w has to be a list
  # to be passed to an apply
  # function
  if (!is.list(w)) 
    w <- list(w)
  ix <- mat2vec.ix(adj, directed, 
                   selfloops)
  # compute the product of all
  # layers to the power of the
  # parameter beta
  product <- apply(sapply(X = 1:length(x), 
                          FUN = function(i) {
                            w[[i]][ix]^x[i]
                          }), MARGIN = 1, FUN = prod)
  # return the value of the right
  # term in eq. 7
  sapply(X = 1:length(x), FUN = function(i) {
    -sum(log(w[[i]][ix]) * xi[ix] * 
           product)/sum(xi[ix] * 
                          product) + sum(adj[ix] * 
                                           log(w[[i]][ix]))/sum(adj[ix])
  })
}


#' Method to predict the expected values of a nrm model
#'
#' @param object nrm object from which to predict
#' @param m integer, the number of edges to be used
#' @param adj optional matrix, the adjacency matrix from which to get the number
#'   of edges
#' @param null optional boolean, is it a null model? default FALSE
#' @param ... other arguments
#' @param multinomial logical. Optional argument. Whether to use multinomial
#'   approximation. If left blank it is selected automatically based on network
#'   size.
#'
#' @return numeric, predicted values from nrm model. (If model is undirected,
#'   only upper.tri of adjacency matrix is returned.)
#' @export
#'
#' @examples
#' data('highschool.predictors')
#' highschool.m <- nrm(w=highschool.predictors[1], adj=contacts.adj, directed=FALSE, selfloops=FALSE)
#' predict(highschool.m, contacts.adj)
#' \donttest{
#' data('highschool.predictors')
#' highschool.m <- nrm(w=highschool.predictors, adj=contacts.adj, directed=FALSE, selfloops=FALSE)
#' predict(highschool.m, contacts.adj)
#' }
predict.nrm <- function(object, 
                        m = NULL, adj = NULL, null = FALSE, 
                        multinomial = NULL, ...) {
  xi <- object$xi
  if (null) {
    omega <- matrix(1, nrow(xi), 
                    ncol(xi))
  } else {
    omega <- matrix(object$omega, 
                    nrow(xi))
  }
  
  selfloops <- object$selfloops
  directed <- object$directed
  ix <- mat2vec.ix(mat = xi, directed = directed, 
                   selfloops = selfloops)
  if (is.null(m) & !is.null(adj)) 
    m <- sum(adj[ix])
  try(if (is.null(m)) 
    stop("Specify number of edges"))
  
  if (is.null(multinomial)) {
    # number of colors and number of
    # links
    if (requireNamespace("BiasedUrn", 
                         quietly = TRUE) && sum(ix) < 
        2000 && (sum(adj[ix])/sum(ix)) < 
        1) {
      multinomial <- FALSE
    } else {
      multinomial <- TRUE
    }
  }
  if (multinomial) {
    return(predict.nrm.multinomial(m, 
                                   xi, omega, ix))
  } else {
    return(predict.nrm.wallenius(m, 
                                 xi, omega, ix))
  }
}

predict.nrm.multinomial <- function(m, 
                                    xi, omega, ix) {
  p <- xi[ix] * omega[ix]/
    sum(xi[ix] * omega[ix])
  return(m * p)
}

predict.nrm.wallenius <- function(m, 
                                  xi, omega, ix) {
  return(BiasedUrn::meanMWNCHypergeo(m = xi[ix], 
                                     n = m, odds = omega[ix]))
}


#' @describeIn as.ghype Map list to ghype
#' @export
#'
as.ghype.nrm <- function(object, ...){
  model <- list(
    'xi'= object$xi,
    'omega' = object$omega,
    'directed' = object$directed,
    'selfloops' = object$selfloops
  )
  class(model) <- 'ghype'
  return(model)
}
