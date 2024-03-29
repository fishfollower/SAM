% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tables.R
\name{faytable}
\alias{faytable}
\alias{faytable.sam}
\title{F-at-age table}
\usage{
faytable(fit, ...)

\method{faytable}{sam}(fit, fleet = which(fit$data$fleetTypes == 0), ...)
}
\arguments{
\item{fit}{a fitted object of class 'sam' as returned from sam.fit}

\item{...}{extra arguments not currently used}

\item{fleet}{the fleet number(s) to return F summed for (default is to return the sum of all residual fleets).}
}
\description{
F-at-age table
}
\details{
...
}
