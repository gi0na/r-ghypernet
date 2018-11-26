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



#' Computes log-likelihood for ghype models.
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
#' @export
logl <- function(adj, xi, omega,
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
        2000 && (sum(adj[ix])) >
        200) {
      multinomial <- FALSE
    } else {
      multinomial <- TRUE
    }
  }
  if (TRUE){ #multinomial) {
    return(logl.multinomial(adj,
                            xi, omega, directed,
                            selfloops))
  } else {
    return(logl.wallenius(adj,
                          xi, omega, directed,
                          selfloops))
  }
}

#' Computes log-likelihood for nrm models with Wallenius Hypergeometric dist.
#' (BiasedUrn)
#'
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param adj  ~~Describe \code{adj} here~~
#' @param xi  ~~Describe \code{xi} here~~
#' @param omega  ~~Describe \code{omega} here~~
#' @param directed  ~~Describe \code{directed} here~~
#' @param selfloops  ~~Describe \code{selfloops} here~~
#' @return scalar
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#'
logl.wallenius <- function(adj,
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


#' Computes approximated log-likelihood for nrm models with multinomial
#' distribution.
#'
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param adj  ~~Describe \code{adj} here~~
#' @param xi  ~~Describe \code{xi} here~~
#' @param omega  ~~Describe \code{omega} here~~
#' @param directed  ~~Describe \code{directed} here~~
#' @param selfloops  ~~Describe \code{selfloops} here~~
#' @return scalar
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#'
logl.multinomial <- function(adj,
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
                   prob = p, log = T)
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


# #' Akaike's An Information Criterion
# #'
# #' Function calculating Akaike's ‘An Information Criterion’
# #' for a ghype fitted model objects for which a log-likelihood
# #' value can be obtained, according to the formula -2log-likelihood + knpar,
# #' where npar represents the number of parameters in the fitted model,
# #' and k = 2 for the usual AIC, or k = log(n) (n being the number of observations)
# #' for the so-called BIC or SBC (Schwarz's Bayesian criterion).
# #'
# #' @param object ghype object
# #' @param ... other parameters passed.
# #' @param k numeric, the penalty per parameter to be used; the default k = 2
# #'
# #' @return If just one object is provided,
# #' a numeric value with the corresponding AIC (or depending on k).
# #' @export
# #'
# AIC.ghype <- function(object, ..., k = 2){
#   2*object$df - 2*object$loglikelihood
# }

# #' Akaike's An Information Criterion
# #'
# #' Function calculating Akaike's ‘An Information Criterion’
# #' for a ghype fitted model objects for which a log-likelihood
# #' value can be obtained, according to the formula -2log-likelihood + knpar,
# #' where npar represents the number of parameters in the fitted model,
# #' and k = 2 for the usual AIC, or k = log(n) (n being the number of observations)
# #' for the so-called BIC or SBC (Schwarz's Bayesian criterion).
# #'
# #' @param object ghype object
# #' @param ... other parameters passed.
# #'
# #' @return If just one object is provided,
# #' a numeric value with the corresponding AIC (or depending on k).
# #'
# #' @export
# #'
# BIC.ghype <- function(object, ...){
#   log(object$m)*object$df - 2*object$loglikelihood
# }
