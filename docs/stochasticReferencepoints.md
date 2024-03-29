# `stochasticReferencepoints`

Estimate stochastic reference points


## Description

The function estimates reference points based on stochastic model forecasts.
 The following reference points are implemented:
 list("\n", "   ", list(list("F=x"), list("F fixed to x, e.g., ", list("\"F=0.3\""), " (NOT IMPLEMENTED YET)")), "\n", "   ", list(list("StatusQuo"), list("F in the last year of the assessment (NOT IMPLEMENTED YET)")), "\n", "   ", list(list("StatusQuo-y"), list("F in the y years before the last in the assessment, e.g., ", list("\"StatusQuo-1\""), " (NOT IMPLEMENTED YET)")), "\n", "   ", list(list("MSY"), list("F that maximizes yield")), "\n", "   ", list(list("0.xMSY"), list("Fs that gives 0.x*100% of MSY, e.g., ", 
    list("\"0.95MSY\""))), "\n", "   ", list(list("Max"), list("F that maximizes yield per recruit")), "\n", "   ", list(list("0.xdYPR"), list("F such that the derivative of yield per recruit is 0.x times the derivative at F=0, e.g., ", list("\"0.1dYPR\""))), "\n", "   ", list(list("0.xSPR"), list("F such that spawners per recruit is 0.x times spawners per recruit at F=0, e.g., ", list("\"0.35SPR\""))), "\n", "   ", list(list("0.xB0"), list("F such that biomass is 0.x times the biomass at F=0, e.g., ", 
    list("\"0.2B0\""))), "\n")


## Usage

```r
stochasticReferencepoints(fit, referencepoints, ...)
list(list("stochasticReferencepoints"), list("sam"))(
  fit,
  referencepoints,
  method = "Q0.5",
  catchType = "catch",
  nYears = 300,
  Frange = c(0, 2),
  nosim = 1000,
  aveYears = max(fit$data$years) + (-9:0),
  selYears = max(fit$data$years),
  newton.control = list(),
  seed = .timeToSeed(),
  formula = ~ibc(F, 5),
  nosim_ci = 200,
  derivedSummarizer = median,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     a sam fit
`referencepoints`     |     a character vector of reference points to estimate (see Details)
`...`     |     additional parameters that can be passed on
`method`     |     estimation method (See Details)
`catchType`     |     catch type: catch, landing, discard
`nYears`     |     Number of years to forecast
`Frange`     |     Range of F values to consider
`nosim`     |     Number of simulations for estimation
`aveYears`     |     Years to average over for biological input
`selYears`     |     Years to average over for selectivity
`newton.control`     |     List of control parameters for optimization
`seed`     |     Seed for simulations
`formula`     |     Formula to estimate optimization criteria as a function of F
`nosim_ci`     |     Number of simulations for bootstrap confidence intervals
`derivedSummarizer`     |     Function to summarize derived per-recruit values


## Details

Reference points can be estimated using these methods:
 
 list("\n", "   ", list(list("Mean"), list("Use least squares to estimate mean equilibrium values")), "\n", "   ", list(list("Q0.x"), list("Use quantile regression to estimate the 0.x quantile of equilibrium values")), "\n", "   ", list(list("Median"), list("Identical to Q0.5")), "\n", "   ", list(list("Mode"), list("(NOT IMPLEMENTED YET)")), "\n") 
 
 To estimate median equilibrium yield, as required by ICES, the method "Q0.5" should be used.


## Value

sam reference point object
 
 reference point object


## Examples

```r
stochasticReferencepoints(fit, c("MSY","0.95MSY","Max","0.35SPR","0.1dYPR"))
```


