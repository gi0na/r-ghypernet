#' Print method for ghype object.
#'
#' @param x ghype model
#' @param ... further arguments passed to or from other methods.
#' @param suppressCall boolean, suppress print of the call
#'
#' @export
#' 
#' @examples
#' data('adj_karate')
#' model <- scm(adj_karate, FALSE, FALSE)
#' print(model)
#'
print.ghype <- function(x, suppressCall = FALSE,
                        ...) {
  # print method for ghype class
  if (!suppressCall) {
    cat("Call:\n")
    print(x$call)
  }

  directed <- 'undirected'
  if(x$directed)
    directed <- 'directed'
  selfloops <- 'no selfloops'
  if(x$selfloops)
    selfloops <- 'selfloops'
  out <- paste('ghype',directed,',',selfloops,'\n')
  cat(out)
  out <- paste(x$n[1], 'vertices,', x$m, 'edges','\n')
  cat(out)
  cat('Loglikelihood:\n')
  cat(x$loglikelihood)
  cat(paste("\ndf:",x$df,'\n'))
}


#' Print method for elements of class \code{'bccm'}.
#'
#' @param x  object of class \code{'bccm'}
#' @param suppressCall  logical, indicating whether to print the call that generated x
#' @param \dots  optional arguments to print or plot methods.
#' @seealso  \code{\link{bccm}}
#' @export
#' @examples 
#' data('adj_karate')
#' data('vertexlabels')
#' bcc.model <- bccm(adj_karate, labels=vertexlabels, directed=FALSE, selfloops=FALSE)
#' print(bcc.model)
#' 
print.bccm <- function(x, suppressCall = FALSE,
                      ...) {
  # print method for ghypeBlock class
  if (!suppressCall) {
    cat("Call:\n")
    print(x$call)
  }

  directed <- 'undirected'
  if(x$directed)
    directed <- 'directed'
  selfloops <- 'no selfloops'
  if(x$selfloops)
    selfloops <- 'selfloops'
  out <- paste('block ghype',directed,',',selfloops,'\n')
  cat(out)
  out <- paste(x$n[1], 'vertices,', x$m, 'edges','\n')
  cat(out)
  cat('Loglikelihood:\n')
  cat(x$loglikelihood)
  cat(paste("\ndf:",x$df,'\n'))

  cat("\nCoefficients:\n")
  cmat <- cbind(x$coef, x$ci[,
                                  3])
  cmat <- cbind(cmat, abs(cmat[,
                               1])/cmat[, 2])
  cmat <- cbind(cmat, 2 * stats::pnorm(-cmat[,
                                             3]))
  colnames(cmat) <- c("Estimate",
                      "Std.Err", "t value", "Pr(>t)")
  stats::printCoefmat(cmat)
}
