# `deterministicReferencepoints`

Function to calculate reference points for the embedded deterministic model of a SAM fit


## Description

The function estimates reference points based on deterministic per-recruit calculations with no process variance.
 The following reference points are implemented:
 list("\n", "   ", list(list("F=x"), list("F fixed to x, e.g., ", list("\"F=0.3\""))), "\n", "   ", list(list("StatusQuo"), list("F in the last year of the assessment")), "\n", "   ", list(list("StatusQuo-y"), list("F in the y years before the last in the assessment, e.g., ", list("\"StatusQuo-1\""))), "\n", "   ", list(list("MSY"), list("F that maximizes yield")), "\n", "   ", list(list("0.xMSY"), list("Fs that gives 0.x*100% of MSY, e.g., ", list("\"0.95MSY\""))), "\n", "   ", list(list("Max"), 
    list("F that maximizes yield per recruit")), "\n", "   ", list(list("0.xdYPR"), list("F such that the derivative of yield per recruit is 0.x times the derivative at F=0, e.g., ", list("\"0.1dYPR\""))), "\n", "   ", list(list("0.xSPR"), list("F such that spawners per recruit is 0.x times spawners per recruit at F=0, e.g., ", list("\"0.35SPR\""))), "\n", "   ", list(list("0.xB0"), list("F such that biomass is 0.x times the biomass at F=0, e.g., ", list("\"0.2B0\""))), "\n")


## Usage

```r
deterministicReferencepoints(fit, referencepoints, ...)
list(list("deterministicReferencepoints"), list("sam"))(
  fit,
  referencepoints,
  catchType = "catch",
  nYears = 300,
  Fsequence = seq(0, 2, len = 50),
  aveYears = max(fit$data$years) + (-9:0),
  selYears = max(fit$data$years),
  biasCorrect = FALSE,
  newton.control = list(),
  run = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     A fitted SAM model
`referencepoints`     |     list of reference points to calculate (See details)
`...`     |     other arguments not used
`catchType`     |     Type of yield to optimize: landing, catch, or discard
`nYears`     |     Number of years in per-recruit calculations
`Fsequence`     |     Sequence of F values for plotting and starting values
`aveYears`     |     Years to average over for biological input
`selYears`     |     Years to average over for selectivity
`biasCorrect`     |     Should bias correction be used in [sdreport](#sdreport) ?
`newton.control`     |     Control arguments passed to the newton optimizer (See [newton](#newton) )
`run`     |     Run estimation? If false, a list of arguments to MakeADFun is returned.


## Value

List of estimated reference points
 
 List of estimated reference points


## Examples

```r
deterministicReferencepoints(fit, c("MSY","0.95MSY","Max","0.35SPR","0.1dYPR","StatusQuo-3"))
```


