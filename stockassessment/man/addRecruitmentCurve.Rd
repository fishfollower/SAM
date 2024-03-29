% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addRecruitmentCurve.R
\name{addRecruitmentCurve}
\alias{addRecruitmentCurve}
\alias{addRecruitmentCurve.sam}
\title{Add stock-recruitment curve to srplot}
\usage{
addRecruitmentCurve(
  fit,
  CI = TRUE,
  col = rgb(0.6, 0, 0),
  cicol = rgb(0.6, 0, 0, 0.3),
  plot = TRUE,
  PI = FALSE,
  picol = rgb(0.6, 0, 0),
  pilty = 2,
  ...
)

\method{addRecruitmentCurve}{sam}(
  fit,
  CI = TRUE,
  col = rgb(0.6, 0, 0),
  cicol = rgb(0.6, 0, 0, 0.3),
  plot = TRUE,
  PI = FALSE,
  picol = rgb(0.6, 0, 0),
  pilty = 2,
  year = NA_real_,
  lastR = NA_real_,
  ...
)
}
\arguments{
\item{fit}{Object to show SR-curve for}

\item{CI}{Add confidence intervals?}

\item{col}{Color of fitted line}

\item{cicol}{Color of confidence intervals}

\item{plot}{Add the curve to a plot?}

\item{PI}{Add prediction intervals?}

\item{picol}{Color of prediction interval line}

\item{pilty}{Line type of prediction interval line}

\item{...}{not used}

\item{year}{Show recruitment calculated conditional on this year (for recruitment functions that depend on year)}

\item{lastR}{Show recruitment calculated conditional on this previous recruitment (for recruitment functions that depend on recruitment the previous year)}
}
\description{
Add stock-recruitment curve to srplot
}
\seealso{
srplot
}
