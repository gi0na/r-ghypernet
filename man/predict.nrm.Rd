% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nrm.R
\name{predict.nrm}
\alias{predict.nrm}
\title{Method to predict the expected values of a nrm model}
\usage{
\method{predict}{nrm}(object, m = NULL, adj = NULL, null = FALSE, multinomial = NULL, ...)
}
\arguments{
\item{object}{nrm object from which to predict}

\item{m}{integer, the number of edges to be used}

\item{adj}{optional matrix, the adjacency matrix from which to get the number
of edges}

\item{null}{optional boolean, is it a null model? default FALSE}

\item{multinomial}{logical. Optional argument. Whether to use multinomial
approximation. If left blank it is selected automatically based on network
size.}

\item{...}{other arguments}
}
\value{
numeric, predicted values from nrm model. (If model is undirected,
  only upper.tri of adjacency matrix is returned.)
}
\description{
Method to predict the expected values of a nrm model
}
\examples{
data('highschool.predictors')
highschool.m <- nrm(w=highschool.predictors[1], adj=contacts.adj, directed=FALSE, selfloops=FALSE)
predict(highschool.m, contacts.adj)
\donttest{
data('highschool.predictors')
highschool.m <- nrm(w=highschool.predictors, adj=contacts.adj, directed=FALSE, selfloops=FALSE)
predict(highschool.m, contacts.adj)
}
}
