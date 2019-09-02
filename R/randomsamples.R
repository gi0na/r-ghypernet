#' Generate random realisations from ghype model.
#'
#' @param nsamples scalar number of realisations
#' @param model ghype model
#' @param m optional scalar, number of edges to draw
#' @param multinomial optional boolean, draw from multinomial?
#' @param seed optional scalar, seed for random sampling.
#' @return
#' list of adjacency matrices.
#'
#' @export
#' 
#' @examples 
#' data('adj_karate')
#' model <- scm(adj_karate, FALSE, FALSE)
#' rghype(1, model)
#'
rghype <- function(nsamples, model, m=NULL, multinomial=NULL, seed=NULL){
  multinomial <- TRUE
  directed <- model$directed
  selfloops <- model$selfloops
  if(is.null(m))
    m <- model$m

  idx <- mat2vec.ix(model$xi,TRUE,TRUE)
  if(length(model$n)==1 | is.null(model$n))
    idx <- mat2vec.ix(model$xi,directed,selfloops)

  n <- model$n
  if(is.null(n))
    n <- nrow(model$xi)

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
  if(all(omega == omega[1]) & m<23331){
    rvec <- t(extraDistr::rmvhyper(nn = nsamples, n = xi, k = m))
  } else{
    if (multinomial) {
      p <- omega*xi/sum(omega*xi)
      rvec <- stats::rmultinom(n = nsamples, size = m, prob = p)
    } else {
      rvec <- cbind(BiasedUrn::rMWNCHypergeo(nran = nsamples, m = xi, n = m, odds = omega))
    }
  }

  graphlist <- lapply(X = 1:ncol(rvec),FUN = function(cls, rvec, directed, selfloops, n){
    vec2mat(vec = rvec[,cls], directed, selfloops,n)},
  rvec=rvec, directed=directed, selfloops=selfloops,n=n)
  if(nsamples==1)
    graphlist <- graphlist[[1]]
  return(graphlist)
}
