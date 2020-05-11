#' Perform AIC forward selection for nrm.
#' 
#' @param adj  the adjacency matrix of the response network
#' @param predictors  list containing the set of predictors as sublists.
#' @param pval  the significance at which computing confidence intervals.
#' @param directed  logical, is the response network directed?
#' @param selfloops  logical, do the response network allows selfloops?
#' @param xi  optional, the possibility matrix \eqn{\Xi}.
#' @param init  optional, initial values passed to the solver to estimate the MLE.
#' @param ncores  optional, number of cores over which parallelise the task.
#' @param \dots  optional arguments to method
#' @return  A nrm object
#' @author  Giona Casiraghi
#' @seealso \code{\link{nrm}}
#' @examples
#' \donttest{
#' data('highschool.predictors')
#' nrm_selection(adj=contacts.adj,predictors=createPredictors(highschool.predictors),
#'   ncores=1,directed=FALSE,selfloops=FALSE)
#'  }
#' @export
nrm_selection <- function(adj, predictors, 
    directed, selfloops, pval = 0.05, 
    xi = NULL, init = NULL, ncores = NULL, 
    ...) UseMethod("nrm_selection", 
    predictors)


#' @describeIn nrm_selection Default method for the nrm stepwise selection.
#' @export
#'
nrm_selection.default <- function(adj, 
                                 predictors, directed, selfloops, 
                                 pval = 0.05, xi = NULL, init = NULL, 
                                 ncores = NULL, ...) {
  stop("Wrong format of predictors: Use createPredictors()")
}


#' @describeIn nrm_selection Method for the nrm stepwise selection when nrmpredictors are passed.
#' @export
#'
nrm_selection.nrmpredictor <- function(adj, 
                                      predictors, directed, selfloops, 
                                      pval = 0.05, xi = NULL, init = NULL, 
                                      ncores = NULL, ...) {
  # Performs a stepwise forward
  # selection with the predictors
  # passed as a list in
  # 'predictors' The default
  # selection is peformed via
  # minimum description length
  # alternatively the
  # selection level can be set to
  # effect size > 0.05
  default.init <- init
  ix <- mat2vec.ix(adj, directed = directed, 
                   selfloops = selfloops)
  M <- sum(adj[ix])
  message("\nEstimating Null model...")
  null.m <- nrm.default(w = list(matrix(1, 
                                        nrow(adj), ncol(adj))), 
                        adj = adj, directed = directed, 
                        selfloops = selfloops, ci = FALSE)
  xi <- null.m$xi
  ww <- predictors
  nms <- c()
  significance <- models <- list()
  DLs <- -null.m$loglikelihood/log(2) + nrow(xi)*(1+directed)/2*
    log(sum(adj[mat2vec.ix(adj,directed,selfloops)]), base = 2)
  mod0 <- null.m
  csR2 <- csR2step <- 0
  message("\nPerforming forward stepwise selection:")
  ## cycle trough all predictors at
  ## each steps choosing the best
  ## one according to AIC
  for (i in 1:length(predictors)) {
    message("\nStep ", i, "...")
    ## select the best predictor
    ## among those in ww and store it
    ## in w
    sel <- nrmChoose(adj = adj, 
                     w.list = ww, xi = xi, 
                     directed = directed, 
                     selfloops = selfloops, 
                     pval = pval, init = init, 
                     ncores = ncores)
    if (is.null(xi)) 
      xi <- sel$xi
    nms <- c(nms, names(predictors)[sel$predictor])
    w <- ww[[sel$predictor]]
    significance <- c(significance, 
                      nr.significance(mod0 = mod0, 
                                      mod1 = sel$model))
    csR2step <- c(csR2step, 
                  coxsnellR2(mod0 = mod0, 
                             mod1 = sel$model, 
                             m = M))
    csR2 <- c(csR2, coxsnellR2(mod0 = null.m, 
                               mod1 = sel$model, m = M))
    DLs <- c(DLs, sel$model$DL)
    mod0 <- sel$model
    ## drop chosen predictor from
    ## predictors
    predictors[[sel$predictor]] <- NULL
    ## create new list of predictors
    ## ww by joining w with the
    ## remaininig predictors in
    ## predictors
    ww <- lapply(X = predictors, 
                 FUN = function(w.new) {
                   c(w, w.new)
                 })
    ## store the intermediate model
    ## computed earlier
    models <- c(models, list(sel$model))
    ## initial values for parameter
    ## estimation in next step
    init <- sel$model$coef
    if (!is.null(default.init)) 
      init <- c(init, default.init)
  }
  ## find best model according to
  ## significance: discard all
  ## predictors with significance
  ## below pval
  message("\nModel estimation concluded.\n")
  AICS <- c(null.m$AIC, sapply(models, 
                               FUN = function(mod) mod$AIC))
  R2S <- c(0, sapply(models, FUN = function(mod) mod$R2))
  selmod <- list(call = match.call(), 
                 M = M, N = nrow(adj), nms = nms, 
                 models = models, significance = c(1, 
                                                   significance), csR2step = csR2step, 
                 mcR2 = R2S, csR2 = csR2, 
                 AIC = AICS, directed = directed, 
                 selfloops = selfloops)
  class(selmod) <- "nrm_selection"
  return(selmod)
}


#' Selects the best set of predictors among the given sets by means of AIC.
#'
#' Computes all the models defined by a list of groups of predictors Returns the
#' best model according to AIC and id of the corresponding predictors in the
#' list The different models are computed in parallel
#'
#' @param adj  adjacency matrix
#' @param w.list  nrmPredictor object. Nested list of predictors to be selected.
#' @param xi  Xi matrix (optional). defaults to scm Xi matrix.
#' @param directed  logical. Is the network directed?
#' @param selfloops  logical. Does the network contain selfloops?
#' @param pval  numeric. the significance at which computing confidence
#'   intervals. defaults to 0.05
#' @param init  initial values for the MLE numerical maximisation. (See
#'   \code{nrm}.)
#' @param ncores  Number of cores for parallelisation of selection process.
#'   (optional) Defaults to number of available cores - 1.
#' @return list containing the best model according to AIC and id of the
#'   corresponding predictors in the list
#'
#' @export
nrmChoose <- function(adj, w.list, 
                      xi = NULL, directed, selfloops, 
                      pval = 0.05, init = NULL, ncores = NULL) {
  # Computes all the models
  # defined by a list of groups of
  # predictors Returns the best
  # model according to AIC and id
  # of the corresponding
  # predictors in the list The
  # different models are computed
  # in parallel
  if (is.null(xi)) {
    xi <- ComputeXi(adj, directed, 
                    selfloops)
  }
  if (is.null(ncores)) {
    ncores <- min(length(w.list), 
                  parallel::detectCores() - 
                    1)
  }
  nr.ms <- parallel::mclapply(X = w.list, 
                              FUN = nrm, adj = adj, xi = xi, 
                              directed = directed, selfloops = selfloops, 
                              pval = pval, significance = FALSE, 
                              init = init, mc.cores = ncores)
  # to.add <- minAIC(nr.ms)
  to.add <- findMDL(nr.ms)
  selected <- list(model = nr.ms[[to.add]], 
                   predictor = to.add, xi = xi)
  class(selected) <- "nrm_selection"
  return(selected)
}


# ' Auxilliary function, finds mininum AIC among different nrm models.
minAIC <- function(nr.ms) {
  # Returns the id of the model
  # with minimal AIC among a
  # selection of models passed as
  # a list.
  min.aic <- Inf
  for (i in 1:length(nr.ms)) {
    if (nr.ms[[i]]$AIC < min.aic) {
      min.aic <- nr.ms[[i]]$AIC
      id <- i
    }
  }
  return(id)
}

# ' Auxilliary function, finds mininum AIC among different nrm models.
# ' 
# '  ~~ A concise (1-5 lines) description of what the function does. ~~
# ' 
# '  ~~ If necessary, more details than the description above ~~
# ' 
# ' @param nr.ms  ~~Describe \code{nr.ms} here~~
# ' @return  
# ' @note  ~~further notes~~
# ' @author  ~~who you are~~
# ' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
# ' @references  ~put references to the literature/web site here ~
# ' @keywords ~kwd1 ~kwd2
# ' @examples
# ' 
# ' 
findMDL <- function(nr.ms) {
  # Returns the id of the model
  # with minimal AIC among a
  # selection of models passed as
  # a list.
  min.dl <- Inf
  for (i in 1:length(nr.ms)) {
    if (nr.ms[[i]]$DL < min.dl) {
      min.dl <- nr.ms[[i]]$DL
      id <- i
    }
  }
  return(id)
}


