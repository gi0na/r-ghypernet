#' Confidence intervals for nrm models.
#' 
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#'  ~~ If necessary, more details than the description above ~~
#' 
#' @param nr.m  ~~Describe \code{nr.m} here~~
#' @param w  ~~Describe \code{w} here~~
#' @param adj  ~~Describe \code{adj} here~~
#' @param pval  ~~Describe \code{pval} here~~
#' @return  
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
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
        pval, "% ci"), "jn")
    ci
}
