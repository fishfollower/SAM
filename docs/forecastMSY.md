# `forecastMSY`

Estimating Fmsy


## Description

Estimating Fmsy


## Usage

```r
forecastMSY(
  fit,
  nYears = 100,
  nlminb.control = list(eval.max = 2000, iter.max = 2000),
  rec.years = c(),
  ave.years = max(fit$data$years) + (-9:0),
  processNoiseF = FALSE,
  ...
)
list(list("forecastMSY"), list("sam"))(
  fit,
  nYears = 100,
  nlminb.control = list(eval.max = 2000, iter.max = 2000, trace = 1),
  rec.years = c(),
  ave.years = max(fit$data$years) + (-9:0),
  processNoiseF = FALSE,
  jacobianHScale = 0.5,
  nCatchAverageYears = 20,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     a SAM fit
`nYears`     |     Number of years to forecast
`nlminb.control`     |     list of control variables for nlminb
`rec.years`     |     Numeric vector of years to use (to calculate mean and standard deviation) for recruitment. An empty vector will use the recruitment model.
`ave.years`     |     vector of years to average for weights, maturity, M and such. Following ICES guidelines, the default is the last 10 years.
`processNoiseF`     |     Should random walk process noise be used for F?
`...`     |     other arguments passed to forecast
`jacobianHScale`     |     Scale step size in jacobian calculation
`nCatchAverageYears`     |     Number of years to average catch over for finding MSY


## Seealso

[forecast](#forecast)  [referencepoints](#referencepoints)


## References

Albertsen, C. M. and Trijoulet, V. (2020) Model-based estimates of reference points in an age-based state-space stock assessment model. Fisheries Research, 230, 105618. doi: 10.1016/j.fishres.2020.105618


