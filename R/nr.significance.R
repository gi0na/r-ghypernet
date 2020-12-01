#' Computes the significance of more complex model against a simpler model by
#' means of a likelihood ratio test.
#' 
#' @param mod0  null nrm model (optional). defaults to the scm model.
#' @param mod1  alternative nrm model, the model to test
#' @param adj  adjacency matrix for which performing the test. (optional) defaults to the matrix used for \code{mod1}.
#' @return  p-value of the lr test mod0 vs mod1
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
