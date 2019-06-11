# `forecast`: forecast function to do shortterm

## Description


 forecast function to do shortterm


## Usage

```r
forecast(fit, fscale = NULL, catchval = NULL, catchval.exact = NULL,
  fval = NULL, nextssb = NULL, landval = NULL, cwF = NULL,
  nosim = 1000, year.base = max(fit$data$years),
  ave.years = max(fit$data$years) + (-4:0),
  rec.years = max(fit$data$years) + (-9:0), label = NULL,
  overwriteSelYears = NULL, deterministic = FALSE,
  processNoiseF = TRUE, customWeights = NULL, customSel = NULL,
  lagR = FALSE, splitLD = FALSE, addTSB = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
```fit```     |     an assessment object of type sam, as returned from the function sam.fit
```fscale```     |     a vector of f-scales. See details.
```catchval```     |     a vector of target catches. See details.
```catchval.exact```     |     a vector of target catches which will be met without noise. See details.
```fval```     |     a vector of target f values. See details.
```nextssb```     |     a vector target SSB values the following year. See details
```landval```     |     a vector of target catches. See details.
```cwF```     |     a vector target custom weighted F values. customWeights must also be specified
```nosim```     |     number of simulations default is 1000
```year.base```     |     starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before
```ave.years```     |     vector of years to average for weights, maturity, M and such
```rec.years```     |     vector of years to use to resample recruitment from
```label```     |     optional label to appear in short table
```overwriteSelYears```     |     if a vector of years is specified, then the average selectivity of those years is used (not recommended)
```deterministic```     |     option to turn all process noise off (not recommended, as it will likely cause bias)
```processNoiseF```     |     option to turn off process noise in F
```customWeights```     |     a vector of same length as number of age groups giving custom weights (currently only used for weighted average of F calculation)
```customSel```     |     supply a custom selection vector that will then be used as fixed selection in all years after the final assessment year (not recommended)
```lagR```     |     if the second youngest age should be reported as recruits
```splitLD```     |     if TRUE the result is split in landing and discards
```addTSB```     |     if TRUE the total stock biomass (TSB) is added

## Details


 There are four ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the last assessment year (base.year), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). The length of the vector specifies how many years forward the scenarios run.


## Value


 an object of type samforecast


