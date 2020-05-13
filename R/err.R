# rmle
err <- function(hatadj, adj){
  sqrt(mean((log(hatadj + 1) - 
               log(adj + 1))^2))
}


#' Computes the Root Mean Squared Error
#'
#' @param model nrm model estimate
#' @param adj original adjacency matrix 
#' @param null logical, whether to compute using null model
#'
#' @return numeric, root mean squared error of residuals of nrm model fit
#' @export
#'
#' @examples
#' 
#' data('highschool.predictors')
#' highschool.m <- nrm(w=highschool.predictors[1], adj=contacts.adj, directed=FALSE, selfloops=FALSE)
#' RMSE(highschool.m, contacts.adj)
#' 

RMSE <- function(model, adj, null = FALSE) {
  return(sqrt(mean(residuals(model, 
                             adj, RMSLE=FALSE, null=null)^2)))
}


#' Computes the Root Mean Squared Logged Error
#'
#' @param model nrm model estimate
#' @param adj original adjacency matrix
#' @param null logical, whether to compute using null model
#'
#' @return numeric, root mean squared logged error of residuals of nrm model fit
#' @export
#'
#' @examples
#' 
#' data('highschool.predictors')
#' highschool.m <- nrm(w=highschool.predictors[1], adj=contacts.adj, directed=FALSE, selfloops=FALSE)
#' RMSLE(highschool.m, contacts.adj)
#' 

RMSLE <- function(model, adj, null = FALSE) {
  return(sqrt(mean(residuals(model, 
                             adj, RMSLE = TRUE, null=null)^2)))
}
