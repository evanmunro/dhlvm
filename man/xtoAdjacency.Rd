% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataAnalysis.R
\name{xtoAdjacency}
\alias{xtoAdjacency}
\title{Converts N x J response matrix to G x L response frequency matrix, when J>1}
\usage{
xtoAdjacency(X, groups)
}
\arguments{
\item{X}{the N x J response matrix}

\item{groups}{N length vector of group membership indictors in 1,..,G}
}
\value{
Y, a G x L response frequency matrix
}
\description{
Converts N x J response matrix to G x L response frequency matrix, when J>1
}
