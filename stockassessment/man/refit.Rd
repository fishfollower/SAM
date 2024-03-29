% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{refit}
\alias{refit}
\title{Re-fit a model from stockassessment.org}
\usage{
refit(fit, newConf, startingValues, ...)
}
\arguments{
\item{fit}{a sam fit or the name of a fit from stockassessment.org}

\item{newConf}{list changes to the configuration}

\item{startingValues}{list of parameter values to use as starting values}

\item{...}{Arguments passed to sam.fit}
}
\value{
A new sam fit
}
\description{
Re-fit a model from stockassessment.org
}
