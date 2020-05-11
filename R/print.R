#' @describeIn ghype
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
  invisible()
}

#' @describeIn bccm
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
  invisible()
}

#' @describeIn nrm
#' Print method for elements of class \code{'nrm'}.
#' 
#' @param x  object of class \code{'nrm'}
#' @param suppressCall  logical, indicating whether to print the call that generated x
#' @param \dots  optional arguments to print or plot methods.
#' @author  Giona Casiraghi
#' @seealso  \code{\link{nrm}}
#' @export
print.nrm <- function(x, suppressCall = FALSE, 
                      ...) {
  # print method for nrm class
  if (!suppressCall) {
    cat("Call:\n")
    print(x$call)
  }
  cat("\nCoefficients:\n")
  cmat <- cbind(x$coef, x$confint[, 
                                  3])
  cmat <- cbind(cmat, abs(cmat[, 
                               1])/cmat[, 2])
  cmat <- cbind(cmat, 2 * stats::pnorm(-cmat[, 
                                             3]))
  colnames(cmat) <- c("Estimate", 
                      "Std.Err", "t value", "Pr(>t)")
  stats::printCoefmat(cmat)
  cat("\nR2:\n")
  print(c(`McFadden R2` = x$R2, 
          `Cox Snell R2` = x$csR2))
  invisible()
}

#' @describeIn nrm_selection
#' Print method for elements of class \code{'nrm_selection'}.
#' 
#' @param x  object of class \code{'nrm_selection'}.
#' @param \dots  optional arguments to print or plot methods.
#' @author  Giona Casiraghi
#' @seealso  \code{nrm_selection}
#' @export
print.nrm_selection <- function(x, 
                                ...) {
  # print method for nrm class
  cat("Call:\n")
  print(x$call)
  id <- which(x$csR2step[-1] < 
                0.05)[1] - 1
  if (id <= 1) 
    id <- 1
  print(x$models[[id]], suppressCall = TRUE)
  aics <- x$AIC[c(1, id, length(x$AIC))]
  es <- c(NA, x$csR2[id], coxsnellR2(x$models[[id]], 
                                     x$models[[length(x$models)]], 
                                     m = x$M))
  names(aics) <- names(es) <- c("null", 
                                "sel", "full")
  cat("\nAIC:\n")
  print(aics[2])
  cat("\nNull AIC and full model AIC:\n")
  out <- cbind(aics[c(1, 3)], 
               es[c(1, 3)])
  colnames(out) <- c("AIC", "effect.s")
  print(out)
  invisible()
}

#' Summary method for elements of class \code{'nrm_selection'}.
#'
#' @param object an object of class 'nrm_selection', usually, a result of a call to \code{nrm_selection}. 
#' @param ... further arguments passed to or from other methods.
#' 
#' @return The function \code{\link{summary.nrm_selection}} computes and
#'   returns a list of summary statistics of the fitted
#'   \code{\link{nrm_selection}} model given in \code{object}.
#' @export
#'
summary.nrm_selection <- function(object, 
                                  ...) {
  # summmary method for nrm class
  results <- cbind(mcR2 = round(object$mcR2, 
                                digits = 4), csR2 = round(object$csR2, 
                                                          digits = 4), AIC = round(object$AIC), 
                   effect.s = round(object$csR2step, 
                                    digits = 4))
  # likelihood ratio tests
  if (length(object$nms) > 0) 
    rownames(results) <- c("-", 
                           object$nms)
  
  ans <- list(object=object,results=results)
  class(ans) <- "summary.nrm_selection"
  ans
}

#' @rdname summary.nrm_selection
#'
#' @param x object of class `summary.nrm_selection` returned by [summary.nrm._selection()].
#' @param ... further arguments passed to or from other methods.
#' @export
print.summary.nrm_selection <- function (x, ...){
  # summmary method for nrm class
  print(x[['object']])
  cat("\n----------------------\n")
  cat("\nAIC selection:\n")
  
  print(x[['results']])
  cat("\nFull model:\n")
  print(x[['object']]$models[[length(x[['object']]$models)]], 
        suppressCall = TRUE)
  
  invisible(x)
}

#' Summary method for elements of class \code{'nrm'}.
#' 
#' Currently it provides the same output as \code{'print.nrm'}
#'
#' @param object an object of class 'nrm', usually, a result of a call to \code{nrm}. 
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function \code{\link{summary.nrm}} computes and
#'   returns a list of summary statistics of the fitted
#'   \code{\link{nrm}} model given in \code{object}.
#'   
#' @export
#'
summary.nrm <- function(object, 
                        ...) {
  # summmary method for nrm class
  ans <- list(object=object)
  class(ans) <- "summary.nrm"
  ans
}

#' @rdname summary.nrm
#'
#' @param x object of class `summary.nrm` returned by [summary.nrm()].
#' @param ... further arguments passed to or from other methods.
#' @export
print.summary.nrm <- function (x, ...){
  print(x[['object']])
  invisible(x)
}

#' Extraction method for coefficients of models of class \code{'nrm'}.
#' 
#' @param object  object of class \code{'nrm'}.
#' @param \dots  optional arguments to print methods.
#' @return coefficients of nrm model.
#' @author  Giona Casiraghi
#' @seealso  \code{\link{nrm}}
#' @export
#' @importFrom stats coef predict residuals
coef.nrm <- function(object, ...) {
  # coef method for nrm class
  object$coef
}
