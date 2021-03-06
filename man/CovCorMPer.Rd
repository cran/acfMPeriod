% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovCorMPer.R
\name{CovCorMPer}
\alias{CovCorMPer}
\title{Robust covariance or correlation matrix from the MPer-ACF}
\usage{
CovCorMPer(x, type = c("correlation", "covariance"))
}
\arguments{
\item{x}{a numeric matrix}

\item{type}{character string giving the type of acf to be computed. Allowed values are "correlation" (the default) or "covariance".}
}
\value{
a numeric matrix
}
\description{
Wrapper that computes the covariance or correlation matrix of \code{x} at lag 0 obtained from the robust MPer-ACF.
}
\examples{
data.set <- cbind(fdeaths, mdeaths)
CovCorMPer(data.set)
}
