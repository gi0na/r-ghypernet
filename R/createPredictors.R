#' Create a nrmpredictor object from passed argument
#'
#' @param predictors the dataframe or list of predictors for to apply nrm model selection
#' @param ... additional parameters passed to the different methods (currently disabled)
#'
#' @return nested list of nrmpredictor class
#' @export
#'
#' @examples
#' data('highschool.predictors')
#' predictors <- createPredictors(highschool.predictors)
#' 
createPredictors <- function(predictors, 
    ...) {
    UseMethod("createPredictors")
}


#' Create a nrmpredictor object from list
#'
#' @param predictors the dataframe or list of predictors for to apply nrm model selection
#' @param ... additional parameters used to creating the predictor object (currently disabled)
#'
#' @return nested list of nrmpredictor class
#' @export
#'
#' @examples
#' data('highschool.predictors')
#' predictors <- createPredictors(highschool.predictors)
#'

createPredictors.list <- function(predictors, 
                                  ...) {
  message("Creating predictors list...")
  predictors <- lapply(X = predictors, 
                       FUN = list)
  for (i in 1:length(predictors)) names(predictors[[i]]) <- names(predictors)[i]
  class(predictors) <- "nrmpredictor"
  return(predictors)
}

# createPredictors.data.frame <-
# function(predictor, ...){
# ##TODO }

# createPredictors.matrix <-
# function(predictor, ...){
# ##TODO }
