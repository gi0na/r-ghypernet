#' Method to compute residuals of nrm models
#'
#' @param object nrm object
#' @param adj odjacency against which to compute residuals
#' @param RMSLE logical, return log residuals? default FALSE
#' @param null logical. use null model?
#' @param ... additional parameters to be passed to and from internal functions.
#'
#' @return numeric vector, residuals of nrm model fit against the original data
#' @export
#'
#' @examples
#' data('highschool.predictors')
#' highschool.m <- nrm(w=highschool.predictors, adj=contacts.adj, directed=FALSE, selfloops=FALSE)
#' residuals(highschool.m, contacts.adj)
#' 
residuals.nrm <- function(object, 
                          adj, RMSLE = FALSE, null = FALSE, 
                          ...) {
  xi <- object$xi
  omega <- object$omega
  selfloops <- object$selfloops
  directed <- object$directed
  
  hatadj <- predict(object = object, 
                    adj = adj, null = null)
  
  ix <- mat2vec.ix(mat = xi, directed = directed, 
                   selfloops = selfloops)
  adj <- adj[ix]
  
  if (!RMSLE) 
    return(adj - hatadj)
  if (RMSLE) 
    return(log(hatadj + 1) - 
             log(adj + 1))
}