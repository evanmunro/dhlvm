% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelAnalysis.R
\name{likelihood}
\alias{likelihood}
\title{Calculate likelihood for a single group of responses}
\usage{
likelihood(data, pi, beta)
}
\arguments{
\item{data}{N xJ response matrix}

\item{pi}{1x K vector of group-specific distribtion over classes}

\item{beta}{J-length list of K x L_j class specific response distributions}
}
\value{
Model ikelihood for provided group of data conditional on provided parameters
}
\description{
Calculate likelihood for a single group of responses
}
