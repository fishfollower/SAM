% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deterministic_referencepoints.R
\name{deterministicReferencepoints}
\alias{deterministicReferencepoints}
\title{Function to calculate reference points for the embedded deterministic model of a SAM fit}
\usage{
deterministicReferencepoints(fit, referencepoints, ...)
}
\arguments{
\item{fit}{A fitted SAM model}

\item{referencepoints}{list of reference points to calculate (See details)}

\item{...}{other arguments not used}
}
\value{
List of estimated reference points
}
\description{
The function estimates reference points based on deterministic per-recruit calculations with no process variance.
The following reference points are implemented:
\describe{
   \item{F=x}{F fixed to x, e.g., \code{"F=0.3"}}
   \item{StatusQuo}{F in the last year of the assessment}
   \item{StatusQuo-y}{F in the y years before the last in the assessment, e.g., \code{"StatusQuo-1"}}
   \item{MSY}{F that maximizes yield}
   \item{0.xMSY}{Fs that gives 0.x*100\% of MSY, e.g., \code{"0.95MSY"}}
   \item{Max}{F that maximizes yield per recruit}
   \item{0.xdYPR}{F such that the derivative of yield per recruit is 0.x times the derivative at F=0, e.g., \code{"0.1dYPR"}}
   \item{0.xSPR}{F such that spawners per recruit is 0.x times spawners per recruit at F=0, e.g., \code{"0.35SPR"}}
   \item{0.xB0}{F such that biomass is 0.x times the biomass at F=0, e.g., \code{"0.2B0"}}
}
}
\examples{
\dontrun{
  deterministicReferencepoints(fit, c("MSY","0.95MSY","Max","0.35SPR","0.1dYPR","StatusQuo-3"))
}
}
