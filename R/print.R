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
}

#' Print method for elements of class \code{'nrm.selection'}.
#' 
#' @param x  object of class \code{'nrm.selection'}.
#' @param \dots  optional arguments to print or plot methods.
#' @return  
#' @author  Giona Casiraghi
#' @seealso  \code{nrm.selection}
#' @export
print.nrm.selection <- function(x, 
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
}

#' Summary method for elements of class \code{'nrm.selection'}.
#'
#' @param object an object of class 'nrm.selection', usually, a result of a call to \code{nrmSelection}. 
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
#' @examples
#'
summary.nrm.selection <- function(object, 
                                  ...) {
  # summmary method for nrm class
  print(object)
  cat("\n----------------------\n")
  cat("\nAIC selection:\n")
  results <- cbind(mcR2 = round(object$mcR2, 
                                digits = 4), csR2 = round(object$csR2, 
                                                          digits = 4), AIC = round(object$AIC), 
                   effect.s = round(object$csR2step, 
                                    digits = 4))
  # likelihood ratio tests
  if (length(object$nms) > 0) 
    rownames(results) <- c("-", 
                           object$nms)
  print(results)
  cat("\nFull model:\n")
  print(object$models[[length(object$models)]], 
        suppressCall = TRUE)
}

#' Summary method for elements of class \code{'nrm'}.
#'
#' @param object an object of class 'nrm', usually, a result of a call to \code{nrm}. 
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
#' @examples
#'
summary.nrm <- function(object, 
                        ...) {
  # summmary method for nrm class
  print(object)
}


#' Print method for coeffiients of models of class \code{'nrm'}.
#' 
#' @param object  object of class \code{'nrm'}.
#' @param \dots  optional arguments to print methods.
#' @return  
#' @author  Giona Casiraghi
#' @seealso  \code{\link{nrm}}
#' @export
#' @importFrom stats coef predict residuals
coef.nrm <- function(object, ...) {
  # coef method for nrm class
  object$coef
}
