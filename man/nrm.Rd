% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nrm.R, R/print.R
\name{nrm}
\alias{nrm}
\alias{nrm.default}
\alias{print.nrm}
\title{Fitting gHypEG regression models for multi-edge networks.}
\usage{
nrm(
  w,
  adj,
  xi = NULL,
  pval = 0.01,
  directed = TRUE,
  selfloops = TRUE,
  regular = FALSE,
  ...
)

\method{nrm}{default}(
  w,
  adj,
  xi = NULL,
  pval = 0.01,
  directed = FALSE,
  selfloops = FALSE,
  regular = FALSE,
  ci = TRUE,
  significance = FALSE,
  null = FALSE,
  init = NULL,
  ...
)

\method{print}{nrm}(x, suppressCall = FALSE, ...)
}
\arguments{
\item{w}{an object of class \code{'list'} containing the predictors layers
(explanatory variables/covariates) of the multiplex, passed as adjacency
matrices. The entries of the list can be named.}

\item{adj}{matrix. The adjacency matrix of the response network (dependent
variable).}

\item{xi}{optional matrix. Passes a non-standard \eqn{\Xi} matrix.}

\item{pval}{the significance level used to compute confidence intervals of
the parameters. Per default, set to 0.01.}

\item{directed}{logical. If \code{TRUE} the response variable is considered
the adjacency matrix of directed graph.  If \code{FALSE} only the upper
triangular of \code{adj} is considered. Default set to FALSE.}

\item{selfloops}{logical. Whether selfloops are allowed. Default set to
FALSE.}

\item{regular}{logical. Whether the gHypEG regression should be performed
with correction of combinatorial effects (\code{TRUE}) or without
(\code{FALSE}).}

\item{\dots}{optional arguments to print or plot methods.}

\item{ci}{logical. Whether to compute confidences for the parameters.
Defaults to \code{TRUE}.}

\item{significance}{logical. Whether to test the model significance against
the null by means of lr-test.}

\item{null}{logical. Is this a null model? Used for internal routines.}

\item{init}{numeric. Vector of initial values used for numerical MLE. If only
a single value is passed, this is repeated to match the number of
predictors in \code{w}.}

\item{x}{object of class \code{'nrm'}}

\item{suppressCall}{logical, indicating whether to print the call that generated x}
}
\value{
nrm returns an object of class 'nrm'.

  The function summary is used to obtain and print a summary and analysis of
  the results. The generic accessory functions coefficients, etc, extract
  various useful features of the value returned by nrm.

  An object of class 'nrm' is a list containing at least the following
  components:

  \item{coef }{a named vector of coefficients.} \item{confint }{a named
  matrix with confidence intervals and standard deviation for each
  coefficient.} \item{omega }{the estimated propensity matrix.} \item{xi
  }{the matrix of possibilities.} \item{loglikelihood }{log-likelihood of the
  estimated model.} \item{AIC }{AIC of the estimated model.} \item{R2 }{Mc
  Fadden pseudo R-squared} \item{csR2 }{Cox and Snells pseudo R-squared}
  \item{significance }{the p-value of the likelihood-ratio test for the
  estimated model against the null.}
}
\description{
nrm is used to fit multi-edge network regression models.
}
\section{Methods (by class)}{
\itemize{
\item \code{nrm(default)}: Default method for nrm

}}
\section{Methods (by generic)}{
\itemize{
\item \code{print(nrm)}: Print method for elements of class \code{'nrm'}.

}}
\examples{

## For a complete example see the vignette

data('highschool.predictors')

highschool.m <- nrm(w=highschool.predictors[1], adj=contacts.adj, directed=FALSE,
  selfloops=FALSE)

highschool.m

\donttest{
data('highschool.predictors')

highschool.m <- nrm(w=highschool.predictors, adj=contacts.adj, directed=FALSE,
  selfloops=FALSE)

highschool.m
}

}
\references{
Casiraghi, Giona. 'Multiplex Network Regression: How do relations
  drive interactions?.' arXiv preprint arXiv:1702.02048 (2017).
}
\seealso{
\code{\link{nrm}}
}
\author{
Giona Casiraghi
}
\keyword{models}
\keyword{multivariate}
\keyword{nonlinear}
\keyword{regression}
\keyword{sna}
