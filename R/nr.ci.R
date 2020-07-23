#' Confidence intervals for nrm models.
#' 
#'  Internal function to compute confidence intervals for estimated parameters of nrm model
#' 
#' @param nr.m  nrm model from which getting coefficients
#' @param w  list of predictors
#' @param adj  adjacency matrix 
#' @param pval  numeric. confidence level
#' @return  matrix reporting values of predictors and confidence bounds
#' 
#' @export
nr.ci <- function(nr.m, w, adj, 
    pval) {
    beta <- nr.m$coef
    jn <- Jn(beta = beta, w = w, 
        xi = nr.m$xi, adj = adj, 
        directed = nr.m$directed, 
        selfloops = nr.m$selfloops)
    jn <- sqrt(diag(solve(jn)))
    # Vectorize(stats::pnorm,
    # vectorize.args =
    # 'lower.tail')(0, mean = beta,
    # sd = jn, lower.tail = beta >
    # 0)
    ci <- cbind(beta - stats::qnorm(pval/2, 
        lower.tail = F) * jn, beta + 
        stats::qnorm(pval/2, lower.tail = F) * 
            jn, jn)
    colnames(ci) <- c(paste(1 - 
        pval, "% ci"), paste(1 - 
        pval, "% ci"), "st_err")
    ci
}
