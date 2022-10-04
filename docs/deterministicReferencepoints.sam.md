# `deterministicReferencepoints.sam`

Function to calculate reference points for the embedded deterministic model of a SAM fit


## Description

Function to calculate reference points for the embedded deterministic model of a SAM fit


## Usage

```r
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
`catchType`     |     landing, catch, or discard
`nYears`     |     Number of years in per-recruit calculations
`Fsequence`     |     Sequence of F values for plotting and starting values
`aveYears`     |     Years to average over for biological input
`selYears`     |     Years to average over for selectivity
`biasCorrect`     |     Should bias correction be used in [sdreport](#sdreport) ?
`newton.control`     |     Control arguments passed to the newton optimizer (See [newton](#newton) )
`run`     |     Run estimation? If false, a list of arguments to MakeADFun is returned.
`...`     |     other arguments not used


## Value

List of estimated reference points


