# `modelforecast`

Model based forecast function


## Description

Model based forecast function
 
 Model based forecast function


## Usage

```r
modelforecast(fit, ...)
list(list("modelforecast"), list("sam"))(
  fit,
  fscale = NULL,
  catchval = NULL,
  fval = NULL,
  nextssb = NULL,
  landval = NULL,
  findMSY = NULL,
  hcr = NULL,
  nosim = NULL,
  year.base = max(fit$data$years),
  ave.years = max(fit$data$years) + (-4:0),
  rec.years = c(),
  label = NULL,
  overwriteSelYears = NULL,
  deterministicF = FALSE,
  processNoiseF = FALSE,
  resampleFirst = !is.null(nosim) && nosim > 0,
  customSel = NULL,
  lagR = FALSE,
  splitLD = FALSE,
  addTSB = FALSE,
  biasCorrect = FALSE,
  returnAllYears = FALSE,
  useUniroot = FALSE,
  nCatchAverageYears = 1,
  returnObj = FALSE,
  hcrConf = numeric(0),
  hcrCurrentSSB = 0,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     an assessment object of type sam, as returned from the function sam.fit
`...`     |     other variables used by the methods
`fscale`     |     a vector of f-scales. See details.
`catchval`     |     a vector of target catches. See details.
`fval`     |     a vector of target f values. See details.
`nextssb`     |     a vector target SSB values the following year. See details
`landval`     |     a vector of target catches. See details.
`findMSY`     |     Should not be used. See [forecastMSY](#forecastmsy) .
`hcr`     |     Should not be used. See [hcr](#hcr) .
`nosim`     |     number of simulations. Not used.
`year.base`     |     starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before
`ave.years`     |     vector of years to average for weights, maturity, M and such
`rec.years`     |     vector of years to use to resample recruitment from. If the vector is empty, the stock recruitment model is used.
`label`     |     optional label to appear in short table
`overwriteSelYears`     |     if a vector of years is specified, then the average selectivity of those years is used (not recommended)
`deterministicF`     |     option to set F variance to (almost) zero (not recommended)
`processNoiseF`     |     option to turn off process noise in F
`resampleFirst`     |     Resample base year when nosim > 0?
`customSel`     |     supply a custom selection vector that will then be used as fixed selection in all years after the final assessment year (not recommended)
`lagR`     |     if the second youngest age should be reported as recruits
`splitLD`     |     if TRUE the result is split in landing and discards
`addTSB`     |     if TRUE the total stock biomass (TSB) is added
`biasCorrect`     |     Do bias correction of reported variables. Can be turned off to reduce running time (not recommended).
`returnAllYears`     |     If TRUE, all years are bias corrected. Otherwise, only forecast years are corrected.
`useUniroot`     |     Use uniroot to find catchval, landval, and nextssb?
`nCatchAverageYears`     |     Should not be used. See [forecastMSY](#forecastmsy) .
`returnObj`     |     Only return TMB object?
`hcrConf`     |     Should not be used. See [hcr](#hcr) .
`hcrCurrentSSB`     |     Should not be used. See [hcr](#hcr) .


## Details

There are four ways to specify a scenario. If e.g. four F values are specified (e.g. fval=c(.1,.2,.3,4)), then the first value is used in the year after the last assessment year (base.year + 1), and the three following in the three following years. Alternatively F's can be specified by a scale, or a target catch. Only one option can be used per year. So for instance to set a catch in the first year and an F-scale in the following one would write catchval=c(10000,NA,NA,NA), fscale=c(NA,1,1,1). If only NA's are specified in a year, the F model is used for forecasting. The length of the vector specifies how many years forward the scenarios run.


## Value

an object of type samforecast


## Seealso

forecast


