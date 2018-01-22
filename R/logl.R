#' Computes log-likelihood for nrm models.
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
#' @param multinomial  ~~Describe \code{multinomial} here~~
#' @return
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
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
        2000 && (sum(adj[ix])/sum(ix)) <
        1) {
      multinomial <- FALSE
    } else {
      multinomial <- TRUE
    }
  }
  if (multinomial) {
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
#' @return
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
#' @return
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



#' Computes log-likelihood ratio for nrm models.
#'
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#'
#'  ~~ If necessary, more details than the description above ~~
#'
#' @param mod0  ~~Describe \code{mod0} here~~
#' @param mod1  ~~Describe \code{mod1} here~~
#' @return val
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' @export
loglratio <- function(mod0, mod1) {
  mod0$loglikelihood - mod1$loglikelihood
}

