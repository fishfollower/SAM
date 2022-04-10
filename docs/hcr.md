# `hcr`

Harvest control rule forecast


## Description

The formula below is used to determine a new F based on the previous SSB.
 
$$F = \left\{\begin{array}{ll}F_{cap} & SSB < B_{cap} \\min\left(Ftarget, \max\left( F_{origin}, (SSB - B_{origin}) \cdot (F_{target} - F_{origin}) / (B_{trigger}-B_{origin}) \right)\right) & SSB \ge B_{origin}\end{array}\right.$$


## Usage

```r
hcr(fit, ...)
hcr(fit, ...)
list(list("hcr"), list("sam"))(
  fit,
  nYears = 20,
  Ftarget,
  Btrigger,
  Forigin = 0,
  Borigin = 0,
  Fcap = 0,
  Bcap = 0,
  nosim = 10000,
  ave.years = max(fit$data$years) + (-4:0),
  rec.years = numeric(0),
  preForecast = list(),
  currentSSB = FALSE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`fit`     |     A SAM fit
`...`     |     additional arguments passed to [modelforecast](#modelforecast)
`nYears`     |     Number of years to forecast
`Ftarget`     |     Target F for high SSB
`Btrigger`     |     SSB that triggers the control rule
`Forigin`     |     F used for SSB = Borigin
`Borigin`     |     Between Blim and Btrigger, F values are selected based on linear interpolation from Forigin to Ftarget
`Fcap`     |     F for SSB < Bcap
`Bcap`     |     SSB for which Fcap is used below
`nosim`     |     Number of simulations to do. If NULL a model forecast based on the Laplace approximation is used
`ave.years`     |     vector of years to average for weights, maturity, M and such
`rec.years`     |     vector of years to use to resample recruitment from. If an empty vector is given, recruitment is based on the fitted model.
`preForecast`     |     list of forecast parameters (i.e., fval, fscale, catchval, landval, or nextssb) to use before the HCR
`currentSSB`     |     if TRUE, SSB at the begining of the control rule year is used. If FALSE, SSB at the begining of the previous year is used.


## Value

hcr object


## Seealso

modelforecast
 
 modelforecast


## Author

Christoffer Moesgaard Albertsen


