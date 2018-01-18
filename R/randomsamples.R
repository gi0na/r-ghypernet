#' Title
#'
#' @param nsamples
#' @param model
#' @param m
#' @param multinomial
#'
#' @return
#' @export
#'
#' @examples
RandomGraph <- function(nsamples, model, m, multinomial=NULL){
  directed <- model$directed
  selfloops <- model$selfloops
  idx <- mat2vec.ix(model$xi,directed,selfloops)
  omega <- model$omega[idx]
  xi <- model$xi[idx]
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
    p <- omega*xi/sum(omega*xi)
    rvec <- rmultinom(n = nsamples, size = m, prob = p)
  } else {
    rvec <- BiasedUrn::rMWNCHypergeo(nran = nsamples, m = xi, n = m, odds = omega)
  }
  return(
    lapply(X = 1:ncol(rvec),FUN = function(cls, rvec, directed, selfloops, n)
    {vec2mat(vec = rvec[,cls], directed, selfloops,n)},
    rvec=rvec, directed=directed, selfloops=selfloops,n=nrow(model$xi)))
}
