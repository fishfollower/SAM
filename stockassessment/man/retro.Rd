% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retro.R
\name{retro}
\alias{retro}
\alias{retro.sam}
\title{retro run}
\usage{
retro(fit, year = NULL, ncores = detectCores(), ...)

\method{retro}{sam}(fit, year = NULL, ncores = detectCores(), ...)
}
\arguments{
\item{fit}{a fitted model object as returned from sam.fit}

\item{year}{either 1) a single integer n in which case runs where all fleets are reduced by 1, 2, ..., n are returned, 2) a vector of years in which case runs where years from and later are excluded for all fleets, and 3 a matrix of years were each column is a fleet and each column corresponds to a run where the years and later are excluded.}

\item{ncores}{the number of cores to attempt to use}

\item{...}{extra arguments to \code{\link{sam.fit}}}
}
\description{
retro run
}
\details{
...
}
