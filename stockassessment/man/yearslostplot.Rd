% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{yearslostplot}
\alias{yearslostplot}
\alias{yearslostplot.default}
\alias{yearslostplot.samforecast}
\alias{yearslostplot.hcr}
\title{SAM years lost to fishing plot}
\usage{
yearslostplot(fit, cause, ...)

\method{yearslostplot}{default}(fit, cause = c("Fishing", "Other", "LifeExpectancy"), ...)

\method{yearslostplot}{samforecast}(fit, cause = c("Fishing", "Other", "LifeExpectancy"), ...)

\method{yearslostplot}{hcr}(fit, cause = c("Fishing", "Other", "LifeExpectancy"), ...)
}
\arguments{
\item{fit}{the object returned from sam.fit}

\item{cause}{Fisning, Other, or LifeExpectancy}

\item{...}{extra arguments transferred to plot including the following: \cr
\code{add} logical, plotting is to be added on existing plot \cr
\code{ci} logical, confidence intervals should be plotted \cr
\code{cicol} color to plot the confidence polygon}
}
\description{
SAM years lost to fishing plot
}
\details{
Plot of years lost to fishing
}
