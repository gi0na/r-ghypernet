#' Computes the significance of more complex model against a simpler model by
#' means of a likelihood ratio test.
#' 
#'  ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#'  ~~ If necessary, more details than the description above ~~
#' 
#' @param mod0  ~~Describe \code{mod0} here~~
#' @param mod1  ~~Describe \code{mod1} here~~
#' @param adj  ~~Describe \code{adj} here~~
#' @return  
#' @note  ~~further notes~~
#' @author  ~~who you are~~
#' @seealso  ~~objects to See Also as \code{\link{help}}, ~~~
#' @references  ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' @export
nr.significance <- function(mod0 = NULL, 
    mod1, adj = NULL) {
    # Perform likelihood-ratio test
    # to quantify significance of
    # more complex model vs simpler
    # model Returns pvalue
    if (is.null(mod0)) {
        mod0 <- nrm.default(w = list(matrix(1, 
            nrow(adj), ncol(adj))), 
            adj = adj, directed = mod1$directed, 
            selfloops = mod1$selfloops, 
            ci = FALSE, significance = FALSE)
        df <- length(mod1$coef)
    } else {
        df <- length(mod1$coef) - 
            length(mod0$coef)
    }
    loglambda <- loglratio(mod0, 
        mod1)
    pval <- stats::pchisq(q = -2 * 
        loglambda, df = df, lower.tail = FALSE)
    return(pval)
}
