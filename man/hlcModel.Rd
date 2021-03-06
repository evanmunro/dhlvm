% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hlcEstimators.R
\name{hlcModel}
\alias{hlcModel}
\title{Run MCMC sampler for LDA-S Model}
\usage{
hlcModel(X, groups, eta, alpha, steps, burn, skip)
}
\arguments{
\item{X}{N x J matrix of responses}

\item{groups}{N x 1 vector of group membership value in 1 to G}

\item{eta}{J length list of K x L_j matrices, Dirichlet prior for beta}

\item{alpha}{G x K matrix, Dirichlet prior for pi}

\item{steps}{Number of MCMC steps}

\item{burn}{Burn-in iterations for sampler}

\item{skip}{Thinning iterations for sampler}
}
\value{
List of full posterior sample for each parameter in the model
}
\description{
Run MCMC sampler for LDA-S Model
}
