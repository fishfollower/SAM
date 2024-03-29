% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotit}
\alias{plotit}
\alias{plotit.sam}
\alias{plotit.samset}
\alias{plotit.samforecast}
\alias{plotit.hcr}
\title{Plot helper}
\usage{
plotit(fit, what, ...)

\method{plotit}{sam}(
  fit,
  what,
  x = fit$data$years,
  ylab = what,
  xlab = "Years",
  ex = numeric(0),
  trans = function(x) x,
  add = FALSE,
  ci = TRUE,
  cicol = gray(0.5, alpha = 0.5),
  addCI = NA,
  drop = 0,
  unnamed.basename = "current",
  xlim = NULL,
  ylim = NULL,
  ylimAdd = NA,
  ...
)

\method{plotit}{samset}(
  fit,
  what,
  x = fit$data$years,
  ylab = what,
  xlab = "Years",
  ex = numeric(0),
  trans = function(x) x,
  add = FALSE,
  ci = TRUE,
  cicol = gray(0.5, alpha = 0.5),
  addCI = rep(FALSE, length(fit)),
  drop = 0,
  unnamed.basename = "current",
  xlim = NULL,
  ...
)

\method{plotit}{samforecast}(
  fit,
  what,
  x = fit$data$years,
  ylab = what,
  xlab = "Years",
  ex = numeric(0),
  trans = function(x) x,
  add = FALSE,
  ci = TRUE,
  cicol = gray(0.5, alpha = 0.5),
  addCI = NA,
  drop = 0,
  unnamed.basename = "current",
  xlim = NULL,
  ylim = NULL,
  ...
)

\method{plotit}{hcr}(
  fit,
  what,
  x = fit$data$years,
  ylab = what,
  xlab = "Years",
  ex = numeric(0),
  trans = function(x) x,
  add = FALSE,
  ci = TRUE,
  cicol = gray(0.5, alpha = 0.5),
  addCI = NA,
  drop = 0,
  unnamed.basename = "current",
  xlim = NULL,
  ...
)
}
\arguments{
\item{fit}{the fitted object from sam.fit of a set of such fits c(fit1,fit2)}

\item{what}{quoted name of object to extract}

\item{...}{extra arguments transferred to plot}

\item{x}{x-values}

\item{ylab}{label on y-axis}

\item{xlab}{label on x-axis}

\item{ex}{extra y's to make room for}

\item{trans}{function to transform values by}

\item{add}{logical, plotting is to be added on existing plot}

\item{ci}{logical, confidence intervals should be plotted}

\item{cicol}{color to plot the confidence polygon}

\item{addCI}{A logical vector indicating if confidence intervals should be plotted for the added fits.}

\item{drop}{number of years to be left unplotted at the end.}

\item{unnamed.basename}{the name to assign an unnamed basefit}

\item{xlim}{xlim for the plot}

\item{ylim}{ylim for the plot}

\item{ylimAdd}{values to add when calculating ylim for the plot}
}
\description{
Plot helper
}
\details{
The basic plotting used bu many of the plotting functions (e.g. ssbplot, fbarplot ...)
}
