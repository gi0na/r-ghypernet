#' Extract Log-Likelihood
#'
#' @param object ghype model.
#' @param ... additional arguments passed to and from internal methods.
#'
#' @return Returns an object of class logLik. This is a number with at least one
#' attribute, "df" (degrees of freedom), giving the number of (estimated) parameters
#' in the model.
#' @export
#'
logLik.ghype <- function(object, ...){
  val <- object$loglikelihood
  attributes(val) <- list(df=object$df, m=object$m, nobs=object$m)
  class(val) <- 'logLik'
  return(val)
}

#' General method to compute log-likelihood for ghype models.
#'
#' @param object  either an adjacency matrix or ghype model If a ghype model is
#'   passed, then `xi`, `omega`, `directed`, `selfloops` are ignored If an
#'   adjacency matrix is passed, then `adj` is ignored
#' @param xi matrix, combinatorial matrix to build ghype model, considered only
#'   if object is an adjacency matrix
#' @param omega matrix, propensity matrix to build ghype model, considered only
#'   if object is an adjacency matrix
#' @param directed boolean, is ghype model directed? considered only if object
#'   is an adjacency matrix
#' @param selfloops boolean, has ghype model selfloops? considered only if
#'   object is an adjacency matrix
#' @param adj optional matrix, adjacency matrix of which to compute
#'   log-likelihood, considered only if object is ghype model If adj is not
#'   passed, and object is a ghype model, the log-likelihood is computed for the
#'   original adjacency matrix stored in object.
#' @param multinomial  optional boolean. Force multinomial approximation? If not
#'   chosen, multinomial chosen for large graphs.
#' @param ... additional parameters passed to and from internal methods
#'   
#' @return
#' loglikelihood value
#'
#' @export
#' 
#' @examples 
#' data('adj_karate')
#' model <- scm(adj_karate, FALSE, FALSE)
#' logl(object = model)
#' new_adj <- adj_karate
#' new_adj[3,4] <- 10
#' logl(object=model, adj=new_adj)
#' 
logl <- function(object, xi=NULL, omega=NULL, directed=NULL, selfloops=NULL, adj = NULL, multinomial = NULL, ...){
  UseMethod('logl')
}

#' @describeIn logl Computes log-likelihood for ghype models from model object
#'
#' @export
#' 
logl.ghype <- function(object, xi=NULL, omega=NULL, directed=NULL, 
                       selfloops=NULL, adj = NULL, multinomial = NULL, ...){
  return(logl_matrix(
              ifelse(test = is.null(adj), yes = list(object$adj), no = list(adj))[[1]],
                     object$xi, object$omega,
                          object$directed, object$selfloops, multinomial = multinomial))
}



#' @describeIn logl Computes log-likelihood for ghype models from adjacency.
#' 
#' @export
#' 
logl.matrix <- function(object, xi=NULL, omega=NULL, directed=NULL, selfloops=NULL, adj = NULL, multinomial = NULL, ...){
  return(logl_matrix(object, xi, omega, directed, selfloops, multinomial = multinomial))
}

# 
# #' Computes log-likelihood for ghype models from adjacency.
# #' 
# #' @param adj  an adjacency matrix
# #' @param xi  a combinatorial matrix
# #' @param omega  a propensity matrix
# #' @param directed  a boolean argument specifying whether object is directed or not.
# #' @param selfloops  a boolean argument specifying whether the model should incorporate selfloops.
# #' @param multinomial  optional boolean. Force multinomial approximation?
# #' If not chosen, multinomial chosen for large graphs.
# #' 
# #' @return
# #' loglikelihood value
# #'
logl_matrix <- function(adj, xi, omega,
                 directed, selfloops, multinomial = NULL) {
  # compute log-likelihood:
  # chooses between Wallenius
  # distribution or multinomial
  # depending on size of the
  # network
  
  # number of colors and number of
  # links
  ix <- mat2vec.ix(adj, directed, selfloops)
  
  ## throw error if all(omega == 0)
  if(all(omega[ix] == 0))
    stop('Not enough pairs with nonzero odds. (all(omega == 0))')
  
  if (is.null(multinomial)) {
    if (requireNamespace("BiasedUrn",
                         quietly = TRUE) && sum(ix) <
        2000 && (sum(adj[ix])) <
        200) {
      multinomial <- FALSE
    } else {
      multinomial <- TRUE
    }
  }

  if( all(omega[ix]==omega[ix][1]) )
    return(logl_hypergeom(adj,
                            xi, directed,
                            selfloops))

  if (multinomial) {
    return(logl_multinomial(adj,
                            xi, omega, directed,
                            selfloops))
  } else {
    return(logl_wallenius(adj,
                          xi, omega, directed,
                          selfloops))
  }
}

# Computes log-likelihood for nrm models with Wallenius Hypergeometric dist.
logl_wallenius <- function(adj,
                           xi, omega, directed, selfloops) {
  # Returns the log-likelihood of
  # the model given by 'xi' and
  # 'omega'.
  if (!requireNamespace("BiasedUrn",
                        quietly = TRUE)) {
    stop("BiasedUrn needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # get indices
  ix <- mat2vec.ix(adj, directed,
                   selfloops)
  # compute likelihood according
  # to wallenius distribution: it
  # uses BiasedUrn package
  log(BiasedUrn::dMWNCHypergeo(x = adj[ix],
                               n = sum(adj[ix]), m = xi[ix],
                               odds = omega[ix]))
}

# Computes log-likelihood for ghype models with Hypergeometric dist.
logl_hypergeom <- function(adj,
                           xi, directed, selfloops) {
  # Returns the log-likelihood of
  # the model given by 'xi' and
  # 'omega'.

  # get indices
  ix <- mat2vec.ix(adj, directed,
                   selfloops)
  # compute likelihood according
  # to hypergeometric distribution
  dmvhyper_base(x = adj[ix],
                       k = sum(adj[ix]), n = xi[ix],
                       log = TRUE)
}


# Computes approximated log-likelihood for ghype models with multinomial
logl_multinomial <- function(adj,
                             xi, omega, directed, selfloops) {
  # Returns the approximated
  # log-likelihood of the model
  # given by 'xi' and 'omega'.
  # get indices
  ix <- mat2vec.ix(adj, directed,
                   selfloops)
  # compute approximated
  # likelihood according to
  # multinomial distribution
  pp <- sum(xi[ix] * omega[ix])
  p <- as.vector(xi[ix] * omega[ix])/pp
  stats::dmultinom(x = as.vector(adj[ix]),
                   prob = p, log = TRUE)
}



#' Compute log-likelihood ratio for ghype models.
#'
#'
#' @param mod0  ghype, null model
#' @param mod1  ghype, alternative model
#' @return scalar, log-likelihood ratio
#' @export
#' 
#' @examples
#' data('adj_karate')
#' sc.model <- scm(adj_karate, FALSE, FALSE)
#' full.model <- ghype(adj_karate, FALSE, FALSE)
#' loglratio(sc.model,full.model)
#' 
loglratio <- function(mod0, mod1) {
  mod0$loglikelihood - mod1$loglikelihood
}

