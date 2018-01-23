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
RandomGraph <- function(nsamples, model, m=NULL, multinomial=NULL, seed=NULL){
  directed <- model$directed
  selfloops <- model$selfloops
  if(is.null(m))
    m <- model$m
  idx <- mat2vec.ix(model$xi,directed,selfloops)
  omega <- model$omega[idx]
  xi <- model$xi[idx]
  if (is.null(multinomial)) {
    # number of colors and number of
    # links
    if (requireNamespace("BiasedUrn",
                         quietly = TRUE) && sum(idx) <
        2000 && (m/sum(idx)) <
        1) {
      multinomial <- FALSE
    } else {
      multinomial <- TRUE
    }
  }
  if(!is.null(seed)){
    set.seed(seed)
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
