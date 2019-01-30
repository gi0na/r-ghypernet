#' Extract Log-Likelihood
#'
#' @param object ghype model.
#' @param ... some methods for this generic function require additional arguments.
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
#' @param object  either an adjacency matrix or ghype model
#' If not chosen, multinomial chosen for large graphs.
#' @param ... extra parameters
#'
#' @return
#' loglikelihood value
#'
#' @export
logl <- function(object, ...){
  UseMethod('logl')
}

#' Computes log-likelihood for ghype models from model object
#'
#' @param object  a ghype model
#' @param adj an optional adjacency matrix for which to compute the loglikelihood
#' @param multinomial  optional boolean. Force multinomial approximation?
#' If not chosen, multinomial chosen for large graphs.
#' @param ... extra parametes
#'
#' @return
#' loglikelihood value
#'
#' @export
logl.ghype <- function(object, adj = NULL, multinomial = NULL, ...){
  return(logl_matrix(
              ifelse(test = is.null(adj), yes = object$adj, no = adj),
                     object$xi, object$omega,
                          object$directed, object$selfloops, multinomial = multinomial))
}



#' Computes log-likelihood for ghype models from adjacency.
#'
#' @param object an adjacency matrix
#' @param xi  a combinatorial matrix
#' @param omega  a propensity matrix
#' @param directed  a boolean argument specifying whether object is directed or not.
#' @param selfloops  a boolean argument specifying whether the model should incorporate selfloops.
#' @param multinomial  optional boolean. Force multinomial approximation?
#' If not chosen, multinomial chosen for large graphs.
#' @param ... extra parameters
#'
#' @return
#' loglikelihood value
#'
#' @export
logl.matrix <- function(object, xi, omega, directed, selfloops, multinomial = NULL, ...){
  return(logl_matrix(object, xi, omega, directed, selfloops, multinomial = multinomial))
}


#' Computes log-likelihood for ghype models from adjacency.
#'
#' @param adj  an adjacency matrix
#' @param xi  a combinatorial matrix
#' @param omega  a propensity matrix
#' @param directed  a boolean argument specifying whether object is directed or not.
#' @param selfloops  a boolean argument specifying whether the model should incorporate selfloops.
#' @param multinomial  optional boolean. Force multinomial approximation?
#' If not chosen, multinomial chosen for large graphs.
#'
#' @return
#' loglikelihood value
#'
logl_matrix <- function(adj, xi, omega,
                 directed, selfloops, multinomial = NULL) {
  # compute log-likelihood:
  # chooses between Wallenius
  # distribution or multinomial
  # depending on size of the
  # network
  if (is.null(multinomial)) {
    # number of colors and number of
    # links
    ix <- mat2vec.ix(adj, directed,
                     selfloops)
    if (requireNamespace("BiasedUrn",
                         quietly = TRUE) && sum(ix) <
        2000 && (sum(adj[ix])) <
        200) {
      multinomial <- FALSE
    } else {
      multinomial <- TRUE
    }
  }

  if( all(omega==omega[1]) )
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
  if (!requireNamespace("extraDistr",
                        quietly = TRUE)) {
    stop("extraDistr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # get indices
  ix <- mat2vec.ix(adj, directed,
                   selfloops)
  # compute likelihood according
  # to wallenius distribution: it
  # uses BiasedUrn package
  extraDistr::dmvhyper(x = adj[ix],
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
loglratio <- function(mod0, mod1) {
  mod0$loglikelihood - mod1$loglikelihood
}

